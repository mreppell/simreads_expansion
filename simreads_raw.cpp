//
// simreads_raw.cpp
//
// Copyright 2012 Darren Kessner
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
//


#include "harp_misc.hpp"
#include "HaplotypeReferenceMulti.hpp"
#include "BAMFile.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <tr1/memory>


using namespace std;
using std::tr1::shared_ptr;
namespace bfs = boost::filesystem;
using namespace boost::lambda; // for _1
namespace ublas = boost::numeric::ublas;
                

struct RawRead
{
    string id;
    string sequence;
    string quality;
};


ostream& operator<<(ostream& os, const RawRead& read)
{
    os << "@" << read.id << endl
       << read.sequence << endl
       << "+" << read.id << endl
       << read.quality << endl;
    return os;
}


class ReadSimulator
{
    public:

    struct Config
    {
        string filename_refseq;
        string filename_refseq_unknown;
        string filename_stem;
        string filename_qh;
        vector<double> haplotype_frequencies;
        double unknown_frequency;
        int unknown_related_haplotype_index; // -1 == none
        double coverage;
        size_t read_length;
        unsigned int seed; 

        Config(const string& filename = "");
        static void print_sample_config();
        void parse_config_file(const string& filename);
        void validate() const;
    }; // Config


    ReadSimulator(const Config& config);

    void simulate_reads() const;


    private:

    Config config_;

    shared_ptr<HaplotypeReferenceMulti> haprefmulti_;
    shared_ptr<HaplotypeReferenceMulti> haprefmulti_unknown_;

    mutable boost::mt19937 gen_;
    boost::random::discrete_distribution<> hapfreq_distribution_;

    vector< boost::random::discrete_distribution<> > base_quality_distributions_;
    void initialize_base_quality_distributions();

    mutable size_t id_count_;
    mutable vector<size_t> haplotype_counts_;
    mutable size_t unknown_count_;

    size_t random_haplotype_index() const {return hapfreq_distribution_(gen_);}
    string add_errors(const string& sequence, const string& quality) const;
    RawRead random_read() const;
};


//
// ReadSimulator::Config
//


ReadSimulator::Config::Config(const string& filename)
:   unknown_frequency(0), unknown_related_haplotype_index(-1),
    coverage(0), read_length(75), 
    seed(static_cast<unsigned int>(std::time(0)))
{
    if (!filename.empty())
        parse_config_file(filename);
}


void ReadSimulator::Config::print_sample_config()
{
    Config config;
    config.filename_refseq = "blah.fasta";
    config.filename_stem = "reads";
    config.haplotype_frequencies.resize(4);
    for (int i=0; i<4; i++) config.haplotype_frequencies[i] = .25; 
    config.coverage = 200;

    ostream& operator<<(ostream& os, const ReadSimulator::Config& config); // forward declaration
    cout << config;
}


void ReadSimulator::Config::parse_config_file(const string& filename)
{
    ifstream is(filename.c_str());
    if (!is) 
        throw runtime_error(("[ReadSimulator::Config::parse_config_file] Unable to open file " + filename).c_str());

    for (string line; getline(is, line);)
    {
        if (line.empty() || line[0] == '#') continue;

        istringstream iss(line);
        string name;
        iss >> name;

        if (name == "filename_refseq")
            iss >> filename_refseq;
        if (name == "filename_refseq_unknown")
            iss >> filename_refseq_unknown;
        else if (name == "filename_stem")
            iss >> filename_stem;
        else if (name == "filename_qh")
            iss >> filename_qh;
        else if (name == "haplotype_frequencies")
            copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(haplotype_frequencies));
        else if (name == "unknown_frequency")
            iss >> unknown_frequency;
        else if (name == "unknown_related_haplotype_index")
            iss >> unknown_related_haplotype_index;
        else if (name == "coverage")
            iss >> coverage;
        else if (name == "seed")
            iss >> seed;
        else if (name == "read_length")
            iss >> read_length;
    }
}


void ReadSimulator::Config::validate() const
{
    if (filename_refseq.empty()) throw runtime_error("[simreads_raw::Config] missing filename_refseq");
    if (filename_stem.empty()) throw runtime_error("[simreads_raw::Config] missing filename_stem");
    if (haplotype_frequencies.empty()) throw runtime_error("[simreads_raw::Config] missing haplotype_frequencies");
    if (coverage == 0) throw runtime_error("[simreads::Config] missing coverage");
}


ostream& operator<<(ostream& os, const ReadSimulator::Config& config)
{
    os << "filename_refseq " << config.filename_refseq << endl;
    os << "filename_refseq_unknown " << config.filename_refseq_unknown << endl;
    os << "filename_stem " << config.filename_stem << endl;
    os << "filename_qh " << config.filename_qh << endl;

    os << "haplotype_frequencies ";
    copy(config.haplotype_frequencies.begin(), config.haplotype_frequencies.end(), ostream_iterator<double>(os, " "));
    os << endl;

    os << "unknown_frequency " << config.unknown_frequency << endl;
    os << "unknown_related_haplotype_index " << config.unknown_related_haplotype_index << endl;
    os << "coverage " << config.coverage << endl;
    os << "seed " << config.seed << endl;

    return os;
}


//
// ReadSimulator
//


ReadSimulator::ReadSimulator(const Config& config)
:   config_(config),
    id_count_(0), unknown_count_(0)
{
    config.validate();

    if (fabs(accumulate(config.haplotype_frequencies.begin(), config.haplotype_frequencies.end(), 0.) - 1) > 1e-6)
        throw runtime_error("[ReadSimulator] Haplotype frequencies do not sum to 1.");

    haprefmulti_ = shared_ptr<HaplotypeReferenceMulti>(new HaplotypeReferenceMulti(config.filename_refseq));
    if (!config.filename_refseq_unknown.empty())
    {
        haprefmulti_unknown_ = shared_ptr<HaplotypeReferenceMulti>(new HaplotypeReferenceMulti(config.filename_refseq_unknown));
        cout << "unknown haplotype count: " << haprefmulti_unknown_-> haplotype_count() << endl;
    }

    if (config.haplotype_frequencies.size() > haprefmulti_->haplotype_count())
        throw runtime_error("[ReadSimulator] Frequency count > haplotype count.");

    haplotype_counts_.resize(config.haplotype_frequencies.size());
    hapfreq_distribution_ = boost::random::discrete_distribution<>(config.haplotype_frequencies);
    initialize_base_quality_distributions();

    gen_.seed(config.seed);
    srand(config.seed);
}


void ReadSimulator::initialize_base_quality_distributions()
{
    if (config_.filename_qh.empty())
    {
        if (config_.read_length != 75) // note: hard-coded
            throw runtime_error("[ReadSimulator::initialize_base_quality_distributions()]] Unsupported read length.");

        extern double mock_02_qual_dist_[41][75];
        
        for (size_t position=0; position<75; ++position)
        {
            vector<double> dist(41);
            for (size_t score=0; score<41; ++score)
                dist[score] = mock_02_qual_dist_[score][position];
		//base_quality_distributions is a vector of distributions, it has 75 entires, one for each position in the read
		//Each distribution is made up of 41 entries, giving the relative probabilities of observing each quality category at that position 
            base_quality_distributions_.push_back(boost::random::discrete_distribution<>(dist));
        }
    }
    else
    {
        ifstream is(config_.filename_qh.c_str());
        if (!is) throw runtime_error(("Unable to open file " + config_.filename_qh).c_str()); 

        cerr << "Reading " << config_.filename_qh << endl;
        vector< vector<double> > qh;

        size_t read_length = 0;

        while (is)
        {
            string buffer;
            getline(is, buffer);
            if (!is) break;
            if (buffer.empty() || buffer[0]=='#') continue;

            qh.push_back(vector<double>());
            istringstream iss(buffer);
            copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(qh.back()));

            if (read_length == 0)
                read_length = qh.back().size();
            else
                if (read_length != qh.back().size())
                    throw runtime_error("qh read length mismatch");
        }

        size_t max_score = qh.size() - 1;

        cerr << "qh read length: " << read_length << endl;
        cerr << "qh max score: " << max_score << endl;

        for (size_t position=0; position<read_length; ++position)
        {
            vector<double> dist(max_score + 1);
            for (size_t score=0; score<=max_score; ++score)
                dist[score] = qh[score][position];
            base_quality_distributions_.push_back(boost::random::discrete_distribution<>(dist));
        }
    }
}


char random_different_base(char notit)
{
    char result = notit;
    while (result == notit)
    {
        double d = rand()/double(RAND_MAX);
        if (d < .25) result = 'A';
        else if (d < .5) result = 'C';
        else if (d < .75) result = 'G';
        else result = 'T';
    }
    return result;
}


string ReadSimulator::add_errors(const string& sequence, const string& quality) const
{
    if (sequence.size() != quality.size()) throw runtime_error("[ReadSimulator::add_errors] Bad vector sizes.");

    string result = sequence;

    for (size_t i=0; i<sequence.size(); ++i)
    {
        double P_error = bampp::illumina_error_probability(quality[i]);

        if (rand()/double(RAND_MAX) < P_error) {
            result[i] = random_different_base(result[i]);
	}
    }
    return result;
}


RawRead ReadSimulator::random_read() const
{
    shared_ptr<HaplotypeReferenceMulti> hapref;
    size_t h = 0;
    string id_tag;

    if (rand()/double(RAND_MAX) < config_.unknown_frequency)
    {
        hapref = haprefmulti_unknown_;
        h = rand() % haprefmulti_unknown_->haplotype_count();
        id_tag = "_unknown_";
        ++unknown_count_;
    }
    else
    {
        hapref = haprefmulti_;
        h = random_haplotype_index();
        id_tag = "_haplotype_";
        ++haplotype_counts_[h];
    }

    size_t haplotype_length = hapref->sequence_length(h);

   // cout << "H-length: " << haplotype_length << " config-length: " << config_.read_length << endl;
    boost::random::uniform_int_distribution<> position_distribution(0, haplotype_length - config_.read_length);
    size_t position = position_distribution(gen_);

    ostringstream id;
    id << "random_read_" << id_count_++ << id_tag << h << "_position_" << position;

    string sequence = hapref->full_sequence(h, position, position + config_.read_length);

    string quality(config_.read_length, '\0');
    for (size_t i=0; i<config_.read_length; ++i)
        quality[i] = base_quality_distributions_[i](gen_); // [0,40]

    RawRead read;
    read.id = id.str();
    read.sequence = add_errors(sequence, quality);

    // standard SAM ASCII encoding (not Illumina)
    read.quality = string(config_.read_length, '\0');
    transform(quality.begin(), quality.end(), read.quality.begin(), _1 + 33);

    return read;
}


// void ReadSimulator::simulate_reads() const
// {
//     bfs::path outdir(".");

//     bfs::ofstream os_seed(outdir / (config_.filename_stem + ".seed"));
//     os_seed << "seed: " << config_.seed << endl;

//     bfs::ofstream os_true_freqs(outdir / (config_.filename_stem + ".true.freqs"));
//     os_true_freqs << "Haplotype Frequencies ";
//     copy(config_.haplotype_frequencies.begin(), config_.haplotype_frequencies.end(), 
//          ostream_iterator<double>(os_true_freqs, " "));
//     os_true_freqs << endl;
//     os_true_freqs.close();

//     shared_ptr<HaplotypeReferenceMulti> temp_hapref;
//     temp_hapref = haprefmulti_;

//     size_t read_count = static_cast<size_t>(config_.coverage * temp_hapref->sequence_length(0) / config_.read_length);
//     os_seed << "read_count: " << read_count << endl; // hack: report this in .seed
//     os_seed.close();

//     bfs::ofstream os_fastq(outdir / (config_.filename_stem + ".reads.fastq"));
//     for (size_t i=0; i<read_count; ++i)
//     {
//         RawRead read = random_read();
//         os_fastq << read;
//     }
//     os_fastq.close();

//     size_t sum_haplotype_counts = accumulate(haplotype_counts_.begin(), haplotype_counts_.end(), 0);
//     if (sum_haplotype_counts + unknown_count_ != id_count_) throw runtime_error("I don't know how to count");

//     bfs::ofstream os_actual_freqs(outdir / (config_.filename_stem + ".actual.freqs"));
//     os_actual_freqs << "Actual Haplotype Frequencies ";
//     for (vector<size_t>::const_iterator it=haplotype_counts_.begin(); it!=haplotype_counts_.end(); ++it)
//         os_actual_freqs << *it/double(sum_haplotype_counts) << " ";
//     os_actual_freqs << endl;
//     os_actual_freqs.close();

// //    bfs::ofstream os_actual_unknown_freq(outdir / (config_.filename_stem + ".actual.unknown_freq"));
// //    os_actual_unknown_freq << double(unknown_count_)/id_count_ << endl;
// //    os_actual_unknown_freq.close();
// }

void ReadSimulator::simulate_reads() const
{
  //    bfs::path outdir(".");

    bfs::ofstream os_seed( config_.filename_stem + ".seed");
    os_seed << "seed: " << config_.seed << endl;

    bfs::ofstream os_true_freqs((config_.filename_stem + ".true.freqs"));
    os_true_freqs << "Haplotype Frequencies ";
    copy(config_.haplotype_frequencies.begin(), config_.haplotype_frequencies.end(), 
         ostream_iterator<double>(os_true_freqs, " "));
    os_true_freqs << endl;
    os_true_freqs.close();

    shared_ptr<HaplotypeReferenceMulti> temp_hapref;
    temp_hapref = haprefmulti_;

    size_t read_count = static_cast<size_t>(config_.coverage);
    //size_t read_count = static_cast<size_t>(config_.coverage * temp_hapref->sequence_length(0) / config_.read_length);
    os_seed << "read_count: " << read_count << endl; // hack: report this in .seed
    os_seed.close();

    bfs::ofstream os_fastq((config_.filename_stem + ".reads.fastq"));
    for (size_t i=0; i<read_count; ++i)
    {
        RawRead read = random_read();
        os_fastq << read;
    }
    os_fastq.close();

    size_t sum_haplotype_counts = accumulate(haplotype_counts_.begin(), haplotype_counts_.end(), 0);
    if (sum_haplotype_counts + unknown_count_ != id_count_) throw runtime_error("I don't know how to count");

    bfs::ofstream os_actual_freqs((config_.filename_stem + ".actual.freqs"));
    os_actual_freqs << "Actual Haplotype Frequencies ";
    for (vector<size_t>::const_iterator it=haplotype_counts_.begin(); it!=haplotype_counts_.end(); ++it)
        os_actual_freqs << *it/double(sum_haplotype_counts) << " ";
    os_actual_freqs << endl;
    os_actual_freqs.close();

//    bfs::ofstream os_actual_unknown_freq(outdir / (config_.filename_stem + ".actual.unknown_freq"));
//    os_actual_unknown_freq << double(unknown_count_)/id_count_ << endl;
//    os_actual_unknown_freq.close();
}


int main(int argc, char* argv[])
{
    try
    {
        if (argc<2)
        {
            cout << "#\n";
            cout << "# Usage: simreads_raw config_filename\n";
            cout << "#\n";
            cout << endl;
            ReadSimulator::Config::print_sample_config();
            return 0;
        }

        const char* filename = argv[1];
        ReadSimulator::Config config(filename);
        cout << config << endl;

        ReadSimulator sim(config);
        sim.simulate_reads();
        return 0;
    }
    catch(exception& e)
    {
        cout << e.what() << endl;
        return 1;
    }
    catch(...)
    {
        cout << "Caught unknown exception.\n";
        return 1;
    }
}


