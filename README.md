# simreads_expansion

This code to expands on Darren Kessner's simreads program (https://bitbucket.org/dkessner/harp) to allow simulation of paired-end reads and paired-end reads from short reference sequences. Simreads was originally created for the Harp software project (Kessner et al (2013) PMID: 23364324), but was expanded for the Karp software project (Reppell and Novembre 2018) to allow simulation of paired-end reads. 

In order to use this code, download the original simreads program, add the .cpp files in this directory to the root directory harp/ where all the other ,cpp files for the program are located. Next, substitue this version of the Jamroot file for that originally included with simreads, and then follow the original installation instructions.

These expansions require the same boost dependencies as the original harp/simreads package.

The new executables in /harp/bin will be:

simreads_raw - Single-end reads with variable lengths 
simreads_raw_pend - Paired-end reads with variable lengths
simreads_sf - Paired-end reads where each forward read begins at the first base of the reference and each paired read begins at the final base of the reference.

Each version takes a configuration file identical to those used by the original simreads program.


