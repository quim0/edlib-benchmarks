/*
 * Copyright (c) 2021 Quim Aguado
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <thread>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <pthread.h>
#include "edlib.h"

#define TIMER_START std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
#define TIMER_STOP std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

#define TIMER_MS std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count()

const char *USAGE_STR = "Usage:\n"
                        "edlib-benchmark <file> <max_seq_len> <num_alignments> "
                        "<threads=1>";

class Sequences {
public:
    size_t seq_len;
    size_t num_alignments;
    char* sequences_buffer;
    int* sequences_len;

    Sequences (char* filepath, int seq_len, int num_alignments) :\
                                                    seq_len(seq_len),
                                                    num_alignments(num_alignments) {
        std::cout << "Sequences object:" << std::endl
                  << "\tFile: " << filepath << std::endl
                  << "\tSequence length: " << seq_len << std::endl
                  << "\tNumber of alignments: " << num_alignments << std::endl;

        std::size_t seq_bytes_to_alloc = ((size_t)num_alignments * (size_t)seq_len * 2L);
        std::cout << "Allocating " << (seq_bytes_to_alloc / (1 << 20))
                  << "MiB of memory to store the sequences" << std::endl;
        try {
            this->sequences_buffer = new char[seq_bytes_to_alloc];
        } catch (std::bad_alloc & exception) {
            std::cerr << "bad_alloc detected: " << exception.what();
            exit(-1);
        }
        memset(this->sequences_buffer, 0, seq_bytes_to_alloc);
        this->sequences_len = new int[(size_t)num_alignments * 2L];

        std::ifstream file(filepath, std::ios::binary | std::ios::ate);
        if (file.fail()) {
            std::cerr << "Could not open file: \"" << filepath << "\"" << std::endl;
            // TODO
            exit(-1);
        }

        std::streamsize size = file.tellg();
        file.seekg(0, std::ios::beg);

        TIMER_START

        std::string line;
        size_t sequences_read = 0;
        while(std::getline(file, line) && sequences_read < (num_alignments*2)) {
            strncpy(this->get_sequence(sequences_read),
                    // +1 to avoid the initial > and <
                    line.c_str() + 1,
                    seq_len);
            this->sequences_len[sequences_read] = line.length() - 1;
            sequences_read++;
        }

        TIMER_STOP
        std::cout << "Read " << sequences_read << " sequences in " << TIMER_MS
                  << "ms." << std::endl;
    };

    ~Sequences () {
        delete [] this->sequences_buffer;
        delete [] this->sequences_len;
    }

    char* get_sequence(size_t n) const {
        // Only for debug purposes
        //if (n >= this->num_alignments*2) {
        //    std::cout << "Trying to read too far... n=" << n << std::endl;
        //    return 0;
        //}
        return this->sequences_buffer + (this->seq_len * n);
    }
};

int worker (const Sequences* seqs, const int tid, const int num_threads) {
    int alignments_to_process = seqs->num_alignments / num_threads;

    if (tid == (num_threads-1)) {
        alignments_to_process += (seqs->num_alignments % num_threads);
    }

    std::cout << "Starting thread " << tid << " to process " <<
              alignments_to_process << " alignments." << std::endl;

    int initial_alignment_idx = (seqs->num_alignments / num_threads) * tid;
    EdlibAlignConfig config;
    config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
    for (int i=0; i<alignments_to_process; i++) {
        EdlibAlignResult result;
        int seq_id_query = i*2 + initial_alignment_idx*2;
        int seq_id_target = seq_id_query + 1;
        char* query = seqs->get_sequence(seq_id_query);
        int query_len = seqs->sequences_len[seq_id_query];
        char* target = seqs->get_sequence(seq_id_target);
        int target_len = seqs->sequences_len[seq_id_target];

        //std::cout << "QUERY(" << query_len << "): " << query << std::endl;
        //std::cout << "TARGET(" << target_len << "): " << target << std::endl;

        result = edlibAlign(query, query_len, target, target_len, config);
        edlibFreeAlignResult(result);
    }

    return 0;
}

int main (int argc, char* argv[]) {
    char* filepath;
    int threads = 1;
    int seq_size = 0;
    int num_alignments;

    if (argc >= 4) {
        filepath = argv[1];
        seq_size = std::atoi(argv[2]);
        num_alignments = std::atoi(argv[3]);
    }
    if (argc == 5) {
        threads = std::atoi(argv[4]);
    }
    if (argc < 4 || argc > 5) {
        std::cerr << USAGE_STR << std::endl;
        return EXIT_FAILURE;
    }

    Sequences sequences(filepath, seq_size, num_alignments);

    std::thread* threads_array = new std::thread[threads-1];

    pthread_t self_thread = pthread_self();
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(0, &cpuset);
    int rc = pthread_setaffinity_np(self_thread, sizeof(cpu_set_t), &cpuset);
    if (rc != 0) {
        std::cerr << "Error setting thread affinity for master\n" << std::endl;
    }

    TIMER_START

    for (int i=0; i<threads-1; i++) {
        threads_array[i] = std::thread(worker, &sequences, i+1, threads);

        // Set thread affinity
        CPU_ZERO(&cpuset);
        CPU_SET((i+1), &cpuset);
        rc = pthread_setaffinity_np(threads_array[i].native_handle(),
                                        sizeof(cpu_set_t), &cpuset);
        if (rc != 0) {
            std::cerr << "Error setting thread affinity for thread " << i
                      << "\n" << std::endl;
        }
    }

    // Master acts as thread 0
    worker(&sequences, 0, threads);

    for (int i=0; i<threads-1; i++) {
        threads_array[i].join();
    }

    TIMER_STOP

    std::cout << num_alignments << " alignments calculated, using " << threads
              << " threads." << std::endl;
    std::cout << "Wall time: " << TIMER_MS << "ms." << std::endl;

    delete [] threads_array;
}
