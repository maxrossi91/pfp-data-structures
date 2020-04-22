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
#include <pfp.hpp>
#include <lce_support.hpp>
#include <sa_support.hpp>

extern "C" {
    #include<gsacak.h>
}

#include <benchmark/benchmark.h>



// std::string test_file = "../data/yeast.fasta";

// static void BM_pfp_sa_queries(benchmark::State &state)
// {
//     // Perform setup here
//     std::string filename = test_file + ".pf.ds";

//     pf_parsing pf;
//     sdsl::load_from_file(pf,filename);

//     pfp_lce_support lce_ds(pf);
//     pfp_sa_support sa_ds(pf);


//     for (auto _ : state)
//     {
//         // This code gets timed
//         for (int i = 0; i < sa_ds.size(); ++i)
//             benchmark::DoNotOptimize(sa_ds.sa(i));
//     }

//     state.counters.insert({{"nops", sa_ds.size() - 1}});
// }
// BENCHMARK(BM_pfp_sa_queries);

// class Setup
// {
//     typedef sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> sdsl_cst_t;
//     typedef pf_parsing pfp_sa_t;

//     Setup()
//     {
//         // call your setup function
//         std::cout << "singleton ctor called only once in the whole program" << std::endl;
//     }

//     static Setup *setup;
//     sdsl_cst_t cst;
//     pf_parsing pf;

//     std::string prev_sdsl = "";
//     std::string prev_pfp = "";

// public:

//     static Setup *PerformSetup()
//     {
//         if (setup == nullptr)
//             setup = new Setup();
//         return setup;
//     }

//     sdsl_cst_t& get_sdsl(std::string s)
//     {
//         if (s.compare(prev_sdsl) != 0)
//         {
//             prev_sdsl = s;
            
//             std::string filename = prev_sdsl + ".sdsl.cst";
//             sdsl::load_from_file(cst, filename);
//         }
//         return cst;
//     }

//     pfp_sa_t& get_pfp(std::string s)
//     {
//         if (s.compare(prev_pfp) != 0)
//         {
//             prev_pfp = s;
            
//             std::string filename = prev_pfp + ".pf.ds";
//             sdsl::load_from_file(pf,filename);
//         }
//         return pf;
//     }

// };

// Setup *Setup::setup = nullptr;

template <class... ExtraArgs>
void BM_sdsl_sa(benchmark::State &st, std::string s)
{

    // Setup *setup = Setup::PerformSetup();
    // auto cst = setup->get_sdsl(s);
    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
    std::string filename = s + ".sdsl.cst";
    sdsl::load_from_file(cst, filename);

    int n_iter = st.range(0);

    if (n_iter == -1 || n_iter > cst.csa.size())
        n_iter = cst.csa.size();
    
    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        for (int i = 0; i < n_iter; ++i)
            benchmark::DoNotOptimize(cst.csa[i]);
    }

}

template <class... ExtraArgs>
void BM_pfp_sa(benchmark::State &st, std::string s)
{
    // Setup *setup = Setup::PerformSetup();
    // auto pf = setup->get_pfp(s);

    pf_parsing<> pf;
    std::string filename = s + ".pf.ds";
    sdsl::load_from_file(pf, filename);

    pfp_sa_support<> pf_sa(pf);

    int n_iter = st.range(0);

    if (n_iter == -1 || n_iter > pf_sa.size())
        n_iter = pf_sa.size();
    
    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        for (int i = 0; i < n_iter; ++i)
            benchmark::DoNotOptimize(pf_sa.sa(i));
    }

}

BENCHMARK_CAPTURE(BM_sdsl_sa, yeast, std::string("../data/yeast.fasta"))->Range(8 << 10, 8 << 20)->Arg(-1);
BENCHMARK_CAPTURE(BM_sdsl_sa, chr19_1, std::string("../../../data/Chr19/chr19.1.fa"))->Range(8 << 10, 8 << 20);
BENCHMARK_CAPTURE(BM_pfp_sa, yeast, std::string("../data/yeast.fasta"))->Range(8 << 10, 8 << 20);
BENCHMARK_CAPTURE(BM_pfp_sa, chr19_1, std::string("../../../data/Chr19/chr19.1.fa"))->Range(8 << 10, 8 << 20);

// Run the benchmark
BENCHMARK_MAIN();
