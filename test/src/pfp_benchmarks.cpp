/* pfp - prefix free parsing lce structure test
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
   \file pfp_benchbarks.cpp
   \brief pfp_benchbarks.cpp build and test prefix-free parsing lce data structure.
   \author Massimiliano Rossi
   \date 19/04/2020
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


inline std::vector<std::pair<uint64_t,uint64_t>> get_queries(size_t interval_size, size_t n_queries)
{
    unsigned seed = 0;
    std::default_random_engine generator(seed);

    std::uniform_int_distribution<uint64_t> distribution(0, interval_size - 1);
    auto dice = std::bind(distribution, generator);

    std::vector<std::pair<uint64_t,uint64_t>> queries;
    for (size_t i = 0; i < n_queries; ++i)
    {
        queries.push_back({dice(),dice()});
    }
    return queries;
}

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
    
    auto queries = get_queries(cst.csa.size(),n_iter);

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        for (auto query: queries)
            benchmark::DoNotOptimize( cst.csa[ query.first ] );
    }

}

template <class... ExtraArgs>
void BM_pfp_sa(benchmark::State &st, std::string s)
{
    // Setup *setup = Setup::PerformSetup();
    // auto pf = setup->get_pfp(s);

    pf_parsing<> pf;
    std::string filename = s + pf.filesuffix();
    sdsl::load_from_file(pf, filename);

    pfp_sa_support<> pf_sa(pf);

    int n_iter = st.range(0);

    if (n_iter == -1 || n_iter > pf_sa.size())
        n_iter = pf_sa.size();

    auto queries = get_queries(pf_sa.size(), n_iter);

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        for (auto query: queries)
            benchmark::DoNotOptimize( pf_sa.sa( query.first ) );
    }

}

template <class... ExtraArgs>
void BM_pfp_sa_custom(benchmark::State &st, std::string s)
{
    // Setup *setup = Setup::PerformSetup();
    // auto pf = setup->get_pfp(s);

    pf_parsing<pfp_wt_custom> pf;
    std::string filename = s + pf.filesuffix();
    sdsl::load_from_file(pf, filename);

    pfp_sa_support<pfp_wt_custom> pf_sa(pf);

    int n_iter = st.range(0);

    if (n_iter == -1 || n_iter > pf_sa.size())
        n_iter = pf_sa.size();

    auto queries = get_queries(pf_sa.size(), n_iter);

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        for (auto query: queries)
            benchmark::DoNotOptimize( pf_sa.sa( query.first ) );
    }

}

// LCP benchmarks
template <class... ExtraArgs>
void BM_sdsl_lcp(benchmark::State &st, std::string s)
{

    // Setup *setup = Setup::PerformSetup();
    // auto cst = setup->get_sdsl(s);
    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
    std::string filename = s + ".sdsl.cst";
    sdsl::load_from_file(cst, filename);

    int n_iter = st.range(0);

    if (n_iter == -1 || n_iter > cst.csa.size())
        n_iter = cst.csa.size();

    auto queries = get_queries(cst.csa.size(), n_iter);

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        for (auto query : queries)
            benchmark::DoNotOptimize(cst.lcp[query.first]);
    }
}

template <class... ExtraArgs>
void BM_pfp_lcp(benchmark::State &st, std::string s)
{
    // Setup *setup = Setup::PerformSetup();
    // auto pf = setup->get_pfp(s);

    pf_parsing<> pf;
    std::string filename = s + pf.filesuffix();
    sdsl::load_from_file(pf, filename);

    pfp_sa_support<> pf_sa(pf);
    pfp_lce_support<> pf_lce(pf);

    int n_iter = st.range(0);

    if (n_iter == -1 || n_iter > pf_sa.size() - 1)
        n_iter = pf_sa.size() - 1;

    auto queries = get_queries(pf_sa.size(), n_iter);

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        for (int i = 0; i < n_iter; ++i)
            benchmark::DoNotOptimize(pf_lce.lce(pf_sa.sa(i), pf_sa.sa(i + 1)));
    }
}

template <class... ExtraArgs>
void BM_pfp_lcp_custom(benchmark::State &st, std::string s)
{
    // Setup *setup = Setup::PerformSetup();
    // auto pf = setup->get_pfp(s);

    pf_parsing<pfp_wt_custom> pf;
    std::string filename = s + pf.filesuffix();
    sdsl::load_from_file(pf, filename);

    pfp_sa_support<pfp_wt_custom> pf_sa(pf);
    pfp_lce_support<pfp_wt_custom> pf_lce(pf);

    int n_iter = st.range(0);

    if (n_iter == -1 || n_iter > pf_sa.size() - 1)
        n_iter = pf_sa.size() - 1;

    auto queries = get_queries(pf_sa.size(), n_iter);

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        for (int i = 0; i < n_iter; ++i)
            benchmark::DoNotOptimize(pf_lce.lce(pf_sa.sa(i), pf_sa.sa(i + 1)));
    }
}

// LCE benchmarks
template <class... ExtraArgs>
void BM_sdsl_lce(benchmark::State &st, std::string s)
{

    // Setup *setup = Setup::PerformSetup();
    // auto cst = setup->get_sdsl(s);
    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
    std::string filename = s + ".sdsl.cst";
    sdsl::load_from_file(cst, filename);

    int n_iter = st.range(0);

    if (n_iter == -1 || n_iter > cst.csa.size() - 1)
        n_iter = cst.csa.size() - 1;

    auto queries = get_queries(cst.csa.size(), n_iter);

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        for (auto query : queries)
            benchmark::DoNotOptimize(cst.depth(cst.lca(cst.select_leaf(cst.csa.isa[query.first] + 1), cst.select_leaf(cst.csa.isa[query.second] + 1))));
    }
}

template <class... ExtraArgs>
void BM_pfp_lce(benchmark::State &st, std::string s)
{
    // Setup *setup = Setup::PerformSetup();
    // auto pf = setup->get_pfp(s);

    pf_parsing<> pf;
    std::string filename = s + pf.filesuffix();
    sdsl::load_from_file(pf, filename);

    pfp_sa_support<> pf_sa(pf);
    pfp_lce_support<> pf_lce(pf);

    int n_iter = st.range(0);

    if (n_iter == -1 || n_iter > pf_sa.size() - 1)
        n_iter = pf_sa.size() - 1;

    auto queries = get_queries(pf_sa.size(), n_iter);

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        for (auto query : queries)
            benchmark::DoNotOptimize(pf_lce.lce(query.first, query.second));
    }
}

template <class... ExtraArgs>
void BM_pfp_lce_custom(benchmark::State &st, std::string s)
{
    // Setup *setup = Setup::PerformSetup();
    // auto pf = setup->get_pfp(s);

    pf_parsing<pfp_wt_custom> pf;
    std::string filename = s + pf.filesuffix();
    sdsl::load_from_file(pf, filename);

    pfp_sa_support<pfp_wt_custom> pf_sa(pf);
    pfp_lce_support<pfp_wt_custom> pf_lce(pf);

    int n_iter = st.range(0);

    if (n_iter == -1 || n_iter > pf_sa.size() - 1)
        n_iter = pf_sa.size() - 1;

    auto queries = get_queries(pf_sa.size(), n_iter);

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        for (auto query : queries)
            benchmark::DoNotOptimize(pf_lce.lce(query.first, query.second));
    }
}

// Thresholds benchmarks
template <class... ExtraArgs>
void BM_sdsl_thresholds(benchmark::State &st, std::string s)
{

    // Setup *setup = Setup::PerformSetup();
    // auto cst = setup->get_sdsl(s);
    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
    std::string filename = s + ".sdsl.cst";
    sdsl::load_from_file(cst, filename);

    int n_iter = cst.csa.size();

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        std::vector<size_t> thresholds;
        std::vector<uint64_t> last_seen(256, 0);
        std::vector<bool> never_seen(256, true);
        size_t i = 0;
        // get the charater in position i
        auto c = (cst.csa.bwt[i] == 0 ? Dollar : cst.csa.bwt[i]);

        while (i < cst.csa.size())
        {

            // If it is not the first run
            if (!never_seen[c])
            {
                // fint the beginning of the previous run of the same character.
                // binary search the threshold
                size_t start = last_seen[c], end = i, mid;
                while (end > start + 1)
                {
                    mid = (start + end) >> 1;
                    // LCP(start,mid)
                    size_t lce_start_mid = cst.depth(cst.lca(cst.select_leaf(start + 1), cst.select_leaf(mid + 1)));
                    // LCP(mid +1, end)
                    size_t lce_mid_end = cst.depth(cst.lca(cst.select_leaf(mid + 1 + 1), cst.select_leaf(end + 1)));
                    if (lce_start_mid > lce_mid_end)
                        start = mid;
                    else
                        end = mid;
                }
                thresholds.push_back(start);
            }

            // Skipp the run
            while (i < cst.csa.size() && c == (cst.csa.bwt[i] == 0 ? Dollar : cst.csa.bwt[i]))
                ++i;

            last_seen[c] = i - 1;
            never_seen[c] = false;
            c = (cst.csa.bwt[i] == 0 ? Dollar : cst.csa.bwt[i]);
        }
    }
}

template <class... ExtraArgs>
void BM_pfp_thresholds(benchmark::State &st, std::string s)
{
    // Setup *setup = Setup::PerformSetup();
    // auto pf = setup->get_pfp(s);

    pf_parsing<> pf;
    std::string filename = s + pf.filesuffix();
    sdsl::load_from_file(pf, filename);

    pfp_sa_support<> sa_ds(pf);
    pfp_lce_support<> lce_ds(pf);

    int n_iter = sa_ds.size();

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        std::vector<size_t> thresholds;
        std::vector<uint64_t> last_seen(256, 0);
        std::vector<bool> never_seen(256, true);
        // get the charater in position i
        size_t i = 0;
        auto sn = (sa_ds.sa(i) + pf.w - 1) % pf.n; // suffix number
        auto p_i = pf.rank_b_p(sn + 1);            // phrase number
        auto id_p_i = pf.pars.p[p_i - 1];          // phrase_id of the phrase that i belongs to.
        size_t occ_in_p_i_in_D = pf.dict.select_b_d(id_p_i) + (sn - pf.select_b_p(p_i));
        auto c = pf.dict.d[occ_in_p_i_in_D];

        while (i < sa_ds.size())
        {
            // If it is not the first run
            if (!never_seen[c])
            {
                // fint the beginning of the next run of the same character.
                // binary search the threshold
                size_t start = last_seen[c], end = i, mid;
                while (end > start + 1)
                {
                    mid = (start + end) >> 1;
                    // LCP(start,mid)
                    size_t lce_start_mid = lce_ds.lce(sa_ds.sa(start), sa_ds.sa(mid));
                    // LCP(mid +1, end)
                    size_t lce_mid_end = lce_ds.lce(sa_ds.sa(mid + 1), sa_ds.sa(end));
                    if (lce_start_mid > lce_mid_end)
                        start = mid;
                    else
                        end = mid;
                }
                thresholds.push_back(start);
            }

            auto next_c = c;
            // Skipp the run
            while (i < sa_ds.size() && c == next_c)
            {
                ++i;
                // get the charater in position i
                auto sn = (sa_ds.sa(i) + pf.w - 1) % pf.n; // suffix number
                auto p_i = pf.rank_b_p(sn + 1);            // phrase number
                auto id_p_i = pf.pars.p[p_i - 1];          // phrase_id of the phrase that s belongs to.
                size_t occ_in_p_i_in_D = pf.dict.select_b_d(id_p_i) + (sn - pf.select_b_p(p_i));
                next_c = pf.dict.d[occ_in_p_i_in_D];
            }
            last_seen[c] = i - 1;
            never_seen[c] = false;
            c = next_c;
        }
    }
}

template <class... ExtraArgs>
void BM_pfp_thresholds_custom(benchmark::State &st, std::string s)
{
    // Setup *setup = Setup::PerformSetup();
    // auto pf = setup->get_pfp(s);

    pf_parsing<pfp_wt_custom> pf;
    std::string filename = s + pf.filesuffix();
    sdsl::load_from_file(pf, filename);

    pfp_sa_support<pfp_wt_custom> sa_ds(pf);
    pfp_lce_support<pfp_wt_custom> lce_ds(pf);

    int n_iter = sa_ds.size();

    st.counters.insert({{"n_iter", n_iter}});
    for (auto _ : st)
    {
        // This code gets timed
        std::vector<size_t> thresholds;
        std::vector<uint64_t> last_seen(256, 0);
        std::vector<bool> never_seen(256, true);
        // get the charater in position i
        size_t i = 0;
        auto sn = (sa_ds.sa(i) + pf.w - 1) % pf.n; // suffix number
        auto p_i = pf.rank_b_p(sn + 1);            // phrase number
        auto id_p_i = pf.pars.p[p_i - 1];          // phrase_id of the phrase that i belongs to.
        size_t occ_in_p_i_in_D = pf.dict.select_b_d(id_p_i) + (sn - pf.select_b_p(p_i));
        auto c = pf.dict.d[occ_in_p_i_in_D];

        while (i < sa_ds.size())
        {
            // If it is not the first run
            if (!never_seen[c])
            {
                // fint the beginning of the next run of the same character.
                // binary search the threshold
                size_t start = last_seen[c], end = i, mid;
                while (end > start + 1)
                {
                    mid = (start + end) >> 1;
                    // LCP(start,mid)
                    size_t lce_start_mid = lce_ds.lce(sa_ds.sa(start), sa_ds.sa(mid));
                    // LCP(mid +1, end)
                    size_t lce_mid_end = lce_ds.lce(sa_ds.sa(mid + 1), sa_ds.sa(end));
                    if (lce_start_mid > lce_mid_end)
                        start = mid;
                    else
                        end = mid;
                }
                thresholds.push_back(start);
            }

            auto next_c = c;
            // Skipp the run
            while (i < sa_ds.size() && c == next_c)
            {
                ++i;
                // get the charater in position i
                auto sn = (sa_ds.sa(i) + pf.w - 1) % pf.n; // suffix number
                auto p_i = pf.rank_b_p(sn + 1);            // phrase number
                auto id_p_i = pf.pars.p[p_i - 1];          // phrase_id of the phrase that s belongs to.
                size_t occ_in_p_i_in_D = pf.dict.select_b_d(id_p_i) + (sn - pf.select_b_p(p_i));
                next_c = pf.dict.d[occ_in_p_i_in_D];
            }
            last_seen[c] = i - 1;
            never_seen[c] = false;
            c = next_c;
        }
    }
}

// LCE benchmarks
BENCHMARK_CAPTURE(BM_sdsl_sa, yeast, std::string("../data/yeast.fasta"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_sa, chr19_1, std::string("../../../data/Chr19/chr19.1.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_sa, chr19_16, std::string("../../../data/Chr19/chr19.16.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_sa, chr19_128, std::string("../../../data/Chr19/chr19.128.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_sa, chr19_512, std::string("../../../data/Chr19/chr19.512.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_sa, salmonella_50, std::string("../../../data/Salmonella/salmonella.50.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_sa, salmonella_500, std::string("../../../data/Salmonella/salmonella.500.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_sa, salmonella_5000, std::string("../../../data/Salmonella/salmonella.5000.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_sa, einstein_en, std::string("../../../data/pizzachili/repcorpus/real/einstein.en.txt"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_sa, world_leaders, std::string("../../../data/pizzachili/repcorpus/real/world_leaders"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_sa, cere, std::string("../../../data/pizzachili/repcorpus/real/cere"))->Arg(10000);

BENCHMARK_CAPTURE(BM_pfp_sa, yeast, std::string("../data/yeast.fasta"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_sa, chr19_1, std::string("../../../data/Chr19/chr19.1.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_sa, chr19_16, std::string("../../../data/Chr19/chr19.16.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_sa, chr19_128, std::string("../../../data/Chr19/chr19.128.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_sa, chr19_512, std::string("../../../data/Chr19/chr19.512.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_sa, salmonella_50, std::string("../../../data/Salmonella/salmonella.50.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_sa, salmonella_500, std::string("../../../data/Salmonella/salmonella.500.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_sa, salmonella_5000, std::string("../../../data/Salmonella/salmonella.5000.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_sa, einstein_en, std::string("../../../data/pizzachili/repcorpus/real/einstein.en.txt"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_sa, world_leaders, std::string("../../../data/pizzachili/repcorpus/real/world_leaders"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_sa, cere, std::string("../../../data/pizzachili/repcorpus/real/cere"))->Arg(10000);

// LCP benchmarks
BENCHMARK_CAPTURE(BM_sdsl_lcp, yeast, std::string("../data/yeast.fasta"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lcp, chr19_1, std::string("../../../data/Chr19/chr19.1.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lcp, chr19_16, std::string("../../../data/Chr19/chr19.16.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lcp, chr19_128, std::string("../../../data/Chr19/chr19.128.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lcp, chr19_512, std::string("../../../data/Chr19/chr19.512.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lcp, salmonella_50, std::string("../../../data/Salmonella/salmonella.50.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lcp, salmonella_500, std::string("../../../data/Salmonella/salmonella.500.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lcp, salmonella_5000, std::string("../../../data/Salmonella/salmonella.5000.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lcp, einstein_en, std::string("../../../data/pizzachili/repcorpus/real/einstein.en.txt"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lcp, world_leaders, std::string("../../../data/pizzachili/repcorpus/real/world_leaders"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lcp, cere, std::string("../../../data/pizzachili/repcorpus/real/cere"))->Arg(10000);

BENCHMARK_CAPTURE(BM_pfp_lcp, yeast, std::string("../data/yeast.fasta"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lcp, chr19_1, std::string("../../../data/Chr19/chr19.1.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lcp, chr19_16, std::string("../../../data/Chr19/chr19.16.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lcp, chr19_128, std::string("../../../data/Chr19/chr19.128.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lcp, chr19_512, std::string("../../../data/Chr19/chr19.512.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lcp, salmonella_50, std::string("../../../data/Salmonella/salmonella.50.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lcp, salmonella_500, std::string("../../../data/Salmonella/salmonella.500.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lcp, salmonella_5000, std::string("../../../data/Salmonella/salmonella.5000.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lcp, einstein_en, std::string("../../../data/pizzachili/repcorpus/real/einstein.en.txt"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lcp, world_leaders, std::string("../../../data/pizzachili/repcorpus/real/world_leaders"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lcp, cere, std::string("../../../data/pizzachili/repcorpus/real/cere"))->Arg(10000);


// LCE benchmarks
BENCHMARK_CAPTURE(BM_sdsl_lce, yeast, std::string("../data/yeast.fasta"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lce, chr19_1, std::string("../../../data/Chr19/chr19.1.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lce, chr19_16, std::string("../../../data/Chr19/chr19.16.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lce, chr19_128, std::string("../../../data/Chr19/chr19.128.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lce, chr19_512, std::string("../../../data/Chr19/chr19.512.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lce, salmonella_50, std::string("../../../data/Salmonella/salmonella.50.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lce, salmonella_500, std::string("../../../data/Salmonella/salmonella.500.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lce, salmonella_5000, std::string("../../../data/Salmonella/salmonella.5000.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lce, einstein_en, std::string("../../../data/pizzachili/repcorpus/real/einstein.en.txt"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lce, world_leaders, std::string("../../../data/pizzachili/repcorpus/real/world_leaders"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_sdsl_lce, cere, std::string("../../../data/pizzachili/repcorpus/real/cere"))->Arg(10000);

BENCHMARK_CAPTURE(BM_pfp_lce, yeast, std::string("../data/yeast.fasta"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lce, chr19_1, std::string("../../../data/Chr19/chr19.1.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lce, chr19_16, std::string("../../../data/Chr19/chr19.16.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lce, chr19_128, std::string("../../../data/Chr19/chr19.128.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lce, chr19_512, std::string("../../../data/Chr19/chr19.512.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lce, salmonella_50, std::string("../../../data/Salmonella/salmonella.50.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lce, salmonella_500, std::string("../../../data/Salmonella/salmonella.500.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lce, salmonella_5000, std::string("../../../data/Salmonella/salmonella.5000.fa"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lce, einstein_en, std::string("../../../data/pizzachili/repcorpus/real/einstein.en.txt"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lce, world_leaders, std::string("../../../data/pizzachili/repcorpus/real/world_leaders"))->Arg(10000);
// BENCHMARK_CAPTURE(BM_pfp_lce, cere, std::string("../../../data/pizzachili/repcorpus/real/cere"))->Arg(10000);

// Thresholds benchmarks
BENCHMARK_CAPTURE(BM_sdsl_thresholds, yeast, std::string("../data/yeast.fasta"));
// BENCHMARK_CAPTURE(BM_sdsl_thresholds, chr19_1, std::string("../../../data/Chr19/chr19.1.fa"));
// BENCHMARK_CAPTURE(BM_sdsl_thresholds, chr19_16, std::string("../../../data/Chr19/chr19.16.fa"));
// BENCHMARK_CAPTURE(BM_sdsl_thresholds, chr19_128, std::string("../../../data/Chr19/chr19.128.fa"));
// BENCHMARK_CAPTURE(BM_sdsl_thresholds, chr19_512, std::string("../../../data/Chr19/chr19.512.fa"));
// BENCHMARK_CAPTURE(BM_sdsl_thresholds, salmonella_50, std::string("../../../data/Salmonella/salmonella.50.fa"));
// BENCHMARK_CAPTURE(BM_sdsl_thresholds, salmonella_500, std::string("../../../data/Salmonella/salmonella.500.fa"));
// BENCHMARK_CAPTURE(BM_sdsl_thresholds, salmonella_5000, std::string("../../../data/Salmonella/salmonella.5000.fa"));
// BENCHMARK_CAPTURE(BM_sdsl_thresholds, einstein_en, std::string("../../../data/pizzachili/repcorpus/real/einstein.en.txt"));
// BENCHMARK_CAPTURE(BM_sdsl_thresholds, world_leaders, std::string("../../../data/pizzachili/repcorpus/real/world_leaders"));
// BENCHMARK_CAPTURE(BM_sdsl_thresholds, cere, std::string("../../../data/pizzachili/repcorpus/real/cere"));

BENCHMARK_CAPTURE(BM_pfp_thresholds, yeast, std::string("../data/yeast.fasta"));
// BENCHMARK_CAPTURE(BM_pfp_thresholds, chr19_1, std::string("../../../data/Chr19/chr19.1.fa"));
// BENCHMARK_CAPTURE(BM_pfp_thresholds, chr19_16, std::string("../../../data/Chr19/chr19.16.fa"));
// BENCHMARK_CAPTURE(BM_pfp_thresholds, chr19_128, std::string("../../../data/Chr19/chr19.128.fa"));
// BENCHMARK_CAPTURE(BM_pfp_thresholds, chr19_512, std::string("../../../data/Chr19/chr19.512.fa"));
// BENCHMARK_CAPTURE(BM_pfp_thresholds, salmonella_50, std::string("../../../data/Salmonella/salmonella.50.fa"));
// BENCHMARK_CAPTURE(BM_pfp_thresholds, salmonella_500, std::string("../../../data/Salmonella/salmonella.500.fa"));
// BENCHMARK_CAPTURE(BM_pfp_thresholds, salmonella_5000, std::string("../../../data/Salmonella/salmonella.5000.fa"));
// BENCHMARK_CAPTURE(BM_pfp_thresholds, einstein_en, std::string("../../../data/pizzachili/repcorpus/real/einstein.en.txt"));
// BENCHMARK_CAPTURE(BM_pfp_thresholds, world_leaders, std::string("../../../data/pizzachili/repcorpus/real/world_leaders"));
// BENCHMARK_CAPTURE(BM_pfp_thresholds, cere, std::string("../../../data/pizzachili/repcorpus/real/cere"));

// Run the benchmark
BENCHMARK_MAIN();
