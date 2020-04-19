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



std::string test_file = "../data/yeast.fasta";

static void BM_pfp_sa_queries(benchmark::State &state)
{
    // Perform setup here
    std::string filename = test_file + ".pf.ds";

    pf_parsing pf;
    sdsl::load_from_file(pf,filename);

    pfp_lce_support lce_ds(pf);
    pfp_sa_support sa_ds(pf);


    for (auto _ : state)
    {
        // This code gets timed
        for (int i = 0; i < sa_ds.size(); ++i)
            benchmark::DoNotOptimize(sa_ds.sa(i));
    }

    state.counters.insert({{"nops", sa_ds.size() - 1}});
}
BENCHMARK(BM_pfp_sa_queries);

static void BM_sdsl_sa_queries(benchmark::State &state)
{
    // Perform setup here
    std::string filename = test_file + ".sdsl.cst";

    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
    sdsl::load_from_file(cst, filename);

    state.counters.insert({{"nops", cst.csa.size() - 1}});
    for (auto _ : state)
    {
        // This code gets timed
        for (int i = 0; i < cst.csa.size(); ++i)
            benchmark::DoNotOptimize(cst.csa[i]);
    }
}
// Register the function as a benchmark
BENCHMARK(BM_sdsl_sa_queries);

// Run the benchmark
BENCHMARK_MAIN();
