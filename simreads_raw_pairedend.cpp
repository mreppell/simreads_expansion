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
      size_t mean_pair_distance;
      size_t max_pair_distance;
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
  boost::random::poisson_distribution<> pair_distance_distribution_;

    vector< boost::random::discrete_distribution<> > base_quality_distributions_;
    void initialize_base_quality_distributions();

    mutable size_t id_count_;
    mutable vector<size_t> haplotype_counts_;
    mutable size_t unknown_count_;

  size_t random_haplotype_index() const {return hapfreq_distribution_(gen_);}
  size_t random_pair_distance() const {return min(config_.max_pair_distance,(size_t) pair_distance_distribution_(gen_));}
  string add_errors(const string& sequence, const string& quality) const;
  pair<RawRead,RawRead> random_read_pair() const;
};


//
// ReadSimulator::Config
//


ReadSimulator::Config::Config(const string& filename)
:   unknown_frequency(0), unknown_related_haplotype_index(-1),
    coverage(0), read_length(75), mean_pair_distance(140), max_pair_distance(300),
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
	else if (name == "mean_pair_distance") 
	  iss >> mean_pair_distance;
	else if (name == "max_pair_distance")
	  iss >> max_pair_distance;
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
    pair_distance_distribution_(config.mean_pair_distance), 
    id_count_(0),
    unknown_count_(0)
    
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
            throw runtime_error("[ReadSimulator::initialize_base_quality_distributions()]] Unsupported read length without qh file.");

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

char opposite_base(char notit) {
  char result;
  if (notit == 'A') 
    result = 'T';
  else if (notit == 'T') 
    result = 'A';
  else if (notit == 'G')
    result = 'C';
  else if (notit == 'C')
    result = 'G';
    else 
      result = notit;

	return result;
}


string ReadSimulator::add_errors(const string& sequence, const string& quality) const
{
    if (sequence.size() != quality.size()) throw runtime_error("[ReadSimulator::add_errors] Bad vector sizes.");

    string result = sequence;

    for (size_t i=0; i<sequence.size(); ++i)
    {
        double P_error = bampp::illumina_error_probability(quality[i]);

        if (rand()/double(RAND_MAX) < P_error)
            result[i] = random_different_base(result[i]);
    }

    return result;
}


pair<RawRead,RawRead> ReadSimulator::random_read_pair() const
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
    
    size_t gap = random_pair_distance();
    boost::random::uniform_int_distribution<> position_distribution(0, haplotype_length - (config_.read_length + gap));
    size_t position1 = position_distribution(gen_);
    size_t position2 = position1 + config_.read_length + gap;

    ostringstream id1;
    ostringstream id2;
    id1 << "random_read_" << id_count_ << id_tag << h << "_position1_" << position1 << "_position2_" << position2;
    id2 << "random_read_" << id_count_++ << id_tag << h << "_position1_" << position1 << "_position2_" << position2;
 
    string sequence1 = hapref->full_sequence(h, position1, position1 + config_.read_length);
    string sequence2 = hapref->full_sequence(h, position2 - config_.read_length,position2);

    reverse(sequence2.begin(),sequence2.end());

    for (size_t i=0; i<sequence2.size(); ++i)
      sequence2[i] = opposite_base(sequence2[i]);
    
    string quality1(config_.read_length, '\0');
    string quality2(config_.read_length, '\0');
    for (size_t i=0; i<config_.read_length; ++i) {
        quality1[i] = base_quality_distributions_[i](gen_); // [0,40]
        quality2[i] = base_quality_distributions_[i](gen_); // [0,40]
    }
    
    RawRead read1;
    read1.id = id1.str();
    read1.sequence = add_errors(sequence1, quality1);

    // standard SAM ASCII encoding (not Illumina)
    read1.quality = string(config_.read_length, '\0');
    transform(quality1.begin(), quality1.end(), read1.quality.begin(), _1 + 33);

    RawRead read2;
    read2.id = id2.str();
    read2.sequence = add_errors(sequence2, quality2);

    // standard SAM ASCII encoding (not Illumina)
    read2.quality = string(config_.read_length, '\0');
    transform(quality2.begin(), quality2.end(), read2.quality.begin(), _1 + 33);

    return make_pair(read1, read2);
}


void ReadSimulator::simulate_reads() const
{
  // bfs::path outdir(config_.filename_stem);

    bfs::ofstream os_seed(config_.filename_stem +  ".seed");
    os_seed << "seed: " << config_.seed << endl;

    bfs::ofstream os_true_freqs(config_.filename_stem  + ".true.freqs");
    os_true_freqs << "Haplotype Frequencies ";
    copy(config_.haplotype_frequencies.begin(), config_.haplotype_frequencies.end(), 
         ostream_iterator<double>(os_true_freqs, " "));
    os_true_freqs << endl;
    os_true_freqs.close();

    shared_ptr<HaplotypeReferenceMulti> temp_hapref;
    temp_hapref = haprefmulti_;

    size_t read_pair_count = static_cast<size_t>(config_.coverage);
    //size_t read_pair_count = static_cast<size_t>(config_.coverage * temp_hapref->sequence_length(0) / (config_.read_length*2));
    os_seed << "read_pair_count: " << read_pair_count << endl; // hack: report this in .seed
    os_seed.close();

    bfs::ofstream os_fastqf(config_.filename_stem + ".reads.forward.fastq");
    bfs::ofstream os_fastqr(config_.filename_stem + ".reads.reverse.fastq");
    for (size_t i=0; i<read_pair_count; ++i)
    {
      pair<RawRead,RawRead> pair_reads = random_read_pair();
      os_fastqf << pair_reads.first;
      os_fastqr << pair_reads.second;
    }
    os_fastqf.close();
    os_fastqr.close();

    size_t sum_haplotype_counts = accumulate(haplotype_counts_.begin(), haplotype_counts_.end(), 0);
    if (sum_haplotype_counts + unknown_count_ != id_count_) throw runtime_error("I don't know how to count");

    bfs::ofstream os_actual_freqs(config_.filename_stem + ".actual.freqs");
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


