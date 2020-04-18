/* pfp - prefix free parsing lce data structure test
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file lce_test.cpp
   \brief lce_test.cpp build and test prefix-free parsing lce data structure.
   \author Massimiliano Rossi
   \date 03/04/2020
*/


#include<iostream>
#include<vector>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>

#include <common.hpp>
#include <strdup.hpp>
#include <pfp.hpp>
#include <lce_support.hpp>
#include <sa_support.hpp>

extern "C" {
    #include<gsacak.h>
}

// // BT-CST
// #include <compressed/PBTRLCSACST.h>

#include <benchmark/benchmark.h>

// std::string test_file = "../data/yeast.fasta";

// class MyFixture : public benchmark::Fixture
// {
// public:
//     sdsl::csa_wt<> csa;
//     sdsl::lcp_wt<> lcp;
//     std::vector<uint32_t> plain_csa;
//     size_t w = 10;
//     pf_parsing pf;
//     pfp_lce_support lce_ds;

//     MyFixture() : pf(test_file, w), lce_ds(pf)
//     {

//         // TEST lce_ds
//         std::vector<char> text;
//         read_fasta_file(test_file.c_str(), text);

//         uint8_t num_bytes = 1;
//         // build cst of the Text
//         sdsl::cache_config cc(false); // do not delete temp files after csa construction
//         sdsl::construct_im(csa, static_cast<const char *>(&text[0]), num_bytes);

//         cc.delete_files = true; // delete temp files after lcp construction
//         sdsl::construct_im(lcp, static_cast<const char *>(&text[0]), num_bytes);

//         plain_csa.resize(csa.size());
//         for(int i = 0; i < csa.size(); ++i){
//             plain_csa[i] = csa[i];
//         }
//         info("Fixture Initialization done");
//     }

// };

// BENCHMARK_F(MyFixture, BM_pfp_lce_queries)
// (benchmark::State &st)
// {
//     for (auto _ : st)
//     {
//         // This code gets timed
//         for (int i = 0; i < csa.size() - 1; ++i)
//             benchmark::DoNotOptimize(lce_ds.lce(plain_csa[i], plain_csa[i + 1]));
//     }
//     info("BM_pfp_lce_queries done");
// }
// BENCHMARK_REGISTER_F(MyFixture, BM_pfp_lce_queries);


std::string test_file = "../data/yeast.fasta";

static void BM_pfp_lcp_queries(benchmark::State &state)
{
    // Perform setup here
    size_t w = 10;
    pf_parsing<> pf(test_file,w);
    pfp_lce_support<> lce_ds(pf);
    pfp_sa_support<> sa_ds(pf);

    verbose("Initialization pfp done");

    // computing size
    size_t size = 0;
    size += sizeof(pf.dict.d[0]) * pf.dict.d.size();
    size += sizeof(pf.dict.saD[0]) * pf.dict.saD.size();
    size += sizeof(pf.dict.isaD[0]) * pf.dict.isaD.size();
    size += sizeof(pf.dict.lcpD[0]) * pf.dict.lcpD.size();
    size += sizeof(pf.dict.daD[0]) * pf.dict.daD.size();
    size += sizeof(pf.dict.colex_daD[0]) * pf.dict.colex_daD.size();
    size += sizeof(pf.pars.p[0]) * pf.pars.p.size();
    size += sizeof(pf.pars.saP[0]) * pf.pars.saP.size();
    size += sizeof(pf.pars.isaP[0]) * pf.pars.isaP.size();
    size += sizeof(pf.pars.lcpP[0]) * pf.pars.lcpP.size();  


    for (auto _ : state)
    {
        // This code gets timed
        for (int i = 0; i < sa_ds.size()-1; ++i)
            benchmark::DoNotOptimize(lce_ds.lce(sa_ds.sa(i), sa_ds.sa(i + 1)));
    }

    state.counters.insert({{"size in B", size}, {"nops", sa_ds.size() - 1}});
    verbose("Queries done");
}
BENCHMARK(BM_pfp_lcp_queries);

static void BM_sdsl_lcp_queries(benchmark::State &state)
{
    // Perform setup here
    std::vector<char> text;
    read_fasta_file(test_file.c_str(), text);

    uint8_t num_bytes = 1;
    // build cst of the Text
    // sdsl::cache_config cc(false); // do not delete temp files after csa construction  
    // sdsl::csa_wt<> csa;
    // sdsl::construct_im(csa, static_cast<const char *>(&text[0]), num_bytes);

    // cc.delete_files = true; // delete temp files after lcp construction
    // sdsl::lcp_support_sada<> lcp;
    // sdsl::construct_im(lcp, static_cast<const char *>(&text[0]), num_bytes);

    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
    sdsl::construct_im(cst, static_cast<const char *>(&text[0]), num_bytes);

    state.counters.insert({{"size in B", size_in_bytes(cst)}, {"nops", cst.csa.size() - 1}});
    verbose("Initialization sdsl done");
    for (auto _ : state)
    {
        // This code gets timed
        for (int i = 1; i < cst.csa.size(); ++i)
            benchmark::DoNotOptimize(cst.lcp[i]);
    }
    verbose("Queries done");
}
// Register the function as a benchmark
BENCHMARK(BM_sdsl_lcp_queries);

static void BM_sdsl_lce_queries(benchmark::State &state)
{
    // Perform setup here
    std::vector<char> text;
    read_fasta_file(test_file.c_str(), text);

    uint8_t num_bytes = 1;

    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
    sdsl::construct_im(cst, static_cast<const char *>(&text[0]), num_bytes);

    state.counters.insert({{"size in B", size_in_bytes(cst)}, {"nops", cst.csa.size() - 1}});
    verbose("Initialization sdsl done");
    for (auto _ : state)
    {
        // This code gets timed
        for (int i = 0; i < cst.csa.size()-1; ++i)
            benchmark::DoNotOptimize(cst.depth(cst.lca(cst.select_leaf(cst.csa.isa[i]+1), cst.select_leaf(cst.csa.isa[i+1]+1))));
    }
    verbose("Queries done");
}
// Register the function as a benchmark
BENCHMARK(BM_sdsl_lce_queries);


// static void BM_bt_cst_lcp_queries(benchmark::State &state)
// {
//     // Perform setup here
//     std::vector<char> text;
//     read_fasta_file(test_file.c_str(), text);

//     uint8_t num_bytes = 1;
//     int r = 2;     //The arity of the BlockTree
//     int mll = 128; // The max length that a BlockTree's leaf could represent

//     PBTRLCSACST *cst = new PBTRLCSACST(text, PBTRLCSACST.PAPER, r, mll, c);

//     state.counters.insert({{"size in B", 0}, {"nops", cst->lcp_rlcsa_.size() - 1}});
//     verbose("Initialization bt-cst done");
//     for (auto _ : state)
//     {
//         // This code gets timed
//         for (int i = 1; i < cst->lcp_rlcsa_.size(); ++i)
//             benchmark::DoNotOptimize(cst->lcp_rlcsa_[i]);
//     }
//     verbose("Queries done");
// }
// // Register the function as a benchmark
// BENCHMARK(BM_sdsl_lcp_queries);

// Run the benchmark
BENCHMARK_MAIN();
