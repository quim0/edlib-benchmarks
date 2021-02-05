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
#include "edlib.h"

#define TIMER_START std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
#define TIMER_STOP std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

#define TIMER_MS std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count()

const char *USAGE_STR = "Usage:\n"
                        "edlib-benchmark <file> <max_seq_len> <num_alignments> "
                        "<threads=1>";

class Sequences {
public:
    int seq_len;
    int num_alignments;
    char* sequences_buffer;
    int* sequences_len;

    Sequences (char* filepath, int seq_len, int num_alignments) :\
                                                    seq_len(seq_len),
                                                    num_alignments(num_alignments) {
        std::cout << "Sequences object:" << std::endl
                  << "\tFile: " << filepath << std::endl
                  << "\tSequence length: " << seq_len << std::endl
                  << "\tNumber of alignments: " << num_alignments << std::endl;

        size_t seq_bytes_to_alloc = (num_alignments * seq_len * 2);
        std::cout << "Allocating " << (seq_bytes_to_alloc / (1<<20))
                  << "MiB of memory to store the sequences" << std::endl;
        try {
            this->sequences_buffer = new char[seq_bytes_to_alloc];
        } catch (std::bad_alloc & exception) {
            std::cerr << "bad_alloc detected: " << exception.what();
            exit(-1);
        }
        memset(this->sequences_buffer, 0, num_alignments * seq_len * 2);
        this->sequences_len = new int[num_alignments * 2];

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
        int sequences_read = 0;
        while(std::getline(file, line) && sequences_read < (num_alignments*2)) {
            strncpy(this->sequences_buffer + (sequences_read * seq_len),
                    // +1 to avoid the initial > and <
                    line.c_str() + 1,
                    seq_len);
            this->sequences_len[sequences_read] = line.length() - 1;
            sequences_read++;
        }

        TIMER_STOP
        std::cout << "Read file in " << TIMER_MS << "ms" << std::endl;
    };

    ~Sequences () {
        delete [] this->sequences_buffer;
    }

    char* get_sequence(int n) const {
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
    for (int i=0; i<alignments_to_process; i++) {
        EdlibAlignResult result;
        EdlibAlignConfig config;
        int seq_id_query = i*2 + initial_alignment_idx;
        int seq_id_target = i*2 + 1 + initial_alignment_idx;
        char* query = seqs->get_sequence(seq_id_query);
        int query_len = seqs->sequences_len[seq_id_query];
        char* target = seqs->get_sequence(seq_id_target);
        int target_len = seqs->sequences_len[seq_id_target];

        //std::cout << "QUERY(" << query_len << "): " << query << std::endl;
        //std::cout << "TARGET(" << target_len << "): " << target << std::endl;

        config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
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

    std::thread* threads_array = new std::thread[threads];

    TIMER_START

    for (int i=0; i<threads; i++) {
        threads_array[i] = std::thread(worker, &sequences, i, threads);
    }

    for (int i=0; i<threads; i++) {
        threads_array[i].join();
    }

    TIMER_STOP

    std::cout << num_alignments << " alignments calculated, using " << threads
              << " threads." << std::endl;
    std::cout << "Wall time: " << TIMER_MS << "ms." << std::endl;

    delete [] threads_array;
}
