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
#include <gtest/gtest.h>
#include <pfp.hpp>
#include <lce_support.hpp>
#include <sa_support.hpp>

extern "C" {
    #include<gsacak.h>
}

// Snippet from https://stackoverflow.com/questions/16491675/how-to-send-custom-message-in-google-c-testing-framework
namespace testing
{
namespace internal
{
enum GTestColor
{
    COLOR_DEFAULT,
    COLOR_RED,
    COLOR_GREEN,
    COLOR_YELLOW
};

extern void ColoredPrintf(GTestColor color, const char *fmt, ...);
} // namespace internal
} // namespace testing
#define PRINTF(...)                                                                        \
    do                                                                                     \
    {                                                                                      \
        testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); \
        testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__);    \
    } while (0)

// C++ stream interface
class TestCout : public std::stringstream
{
public:
    ~TestCout()
    {
        PRINTF("%s", str().c_str());
    }
};

#define TEST_COUT TestCout()

//*************************************************************************************

std::string test_file;

TEST(lce_construct_test, paper_example)
{

    std::vector<char> text = {'G','A','T','T','A','C','A','T','#',
                              'G','A','T','A','C','A','T','#',
                              'G','A','T','T','A','G','A','T','A','#','#'};
    std::vector<std::string> dict{"#GATTAC", "ACAT#", "AGATA##", "T#GATAC", "T#GATTAG"};
    std::vector<uint32_t> parse{1, 2, 4, 2, 5, 3, 0};
    std::vector<uint32_t> indices{0, 1, 2, 3, 4};
    std::vector<uint8_t> dict2 = {'#', '#', 'G', 'A', 'T', 'T', 'A', 'C', EndOfWord,
                                  'A', 'C', 'A', 'T', '#', EndOfWord,
                                  'A', 'G', 'A', 'T', 'A', '#', '#', EndOfWord,
                                  'T', '#', 'G', 'A', 'T', 'A', 'C', EndOfWord,
                                  'T', '#', 'G', 'A', 'T', 'T', 'A', 'G', EndOfWord, EndOfDict};
    sdsl::bit_vector b_d = {1, 0, 0, 0, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    const size_t n_phrases = 5;
    const std::vector<size_t> phrase_length{8, 5, 7, 7, 8};
    sdsl::bit_vector b_p = {1, 0, 0, 0, 0, 0,  
                            1, 0, 0,
                            1, 0, 0, 0, 0, 
                            1, 0, 0, 
                            1, 0, 0, 0, 0, 0, 
                            1, 0, 0, 0, 0};

    std::vector<uint32_t> frequencies{0, 1, 2, 1, 1, 1};
    size_t w = 2;



    TEST_COUT << "Begin paper test" << std::endl;


    pf_parsing pf(dict2,parse,frequencies, w);
    TEST_COUT << "Pfp built" << std::endl;

    // TEST n
    EXPECT_EQ(pf.n, text.size());
    TEST_COUT << "Test n" << std::endl;

    // TEST n_phrases
    EXPECT_EQ(pf.dict.n_phrases(), n_phrases);
    TEST_COUT << "Test n_phrases" << std::endl;

    // TEST b_d
    for(size_t i = 0; i < pf.dict.b_d.size(); ++i){
        EXPECT_EQ(pf.dict.b_d[i], b_d[i]) << "at position: " << i;
    }
    TEST_COUT << "Test b_d" << std::endl;

    // TEST phrase_length
    for (size_t i = 0; i < phrase_length.size(); ++i)
    {
        EXPECT_EQ(pf.dict.length_of_phrase(i+1), phrase_length[i]);
    }
    TEST_COUT << "Test phrase_length" << std::endl;

    // TEST b_p
    for (size_t i = 0; i < pf.b_p.size(); ++i)
    {
        EXPECT_EQ(pf.b_p[i], b_p[i]) << "at position: " << i;
    }
    TEST_COUT << "Test b_p" << std::endl;

    pfp_lce_support lce_ds(pf);
    pfp_sa_support sa_ds(pf);

    // TEST lce_ds
    std::vector<char> tmp_text = {'#', 'G', 'A', 'T', 'T', 'A', 'C', 'A', 'T', '#',
                                  'G', 'A', 'T', 'A', 'C', 'A', 'T', '#',
                                  'G', 'A', 'T', 'T', 'A', 'G', 'A', 'T', 'A', 0};
    // build lcp of the Text
    // uint8_t num_bytes = 1;
    // build cst of the Text
    // verbose("Computing CSA of the text");
    // sdsl::cache_config cc(false); // do not delete temp files after csa construction
    // sdsl::csa_wt<> csa;
    // sdsl::construct_im(csa, static_cast<const char *>(&tmp_text[0]), num_bytes);

    // verbose("Computing LCP of the text");
    // cc.delete_files = true; // delete temp files after lcp construction
    // sdsl::lcp_wt<> lcp;
    // sdsl::construct_im(lcp, static_cast<const char *>(&tmp_text[0]), num_bytes);

    // verbose("Testing LCE ds");
    // for (int i = 1; i < text.size() - 1; ++i)
    // {
    //     EXPECT_EQ(lce_ds.lce(sa_ds.sa(i), sa_ds.sa(i + 1)), lcp[i + 1]) << "At positions: " << sa_ds.sa(i) << " " << sa_ds.sa(i + 1);
    // }

    uint8_t num_bytes = 1;
    // build cst of the Text
    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
    sdsl::construct_im(cst, static_cast<const char *>(&text[0]), num_bytes);

    TEST_COUT << "Testing SA ISA ds" << std::endl;
    for (int i = 0; i < text.size() - 1; ++i)
    {
        EXPECT_EQ(cst.csa[cst.csa.isa[i]], i) << "At position: " << i;
        EXPECT_EQ(cst.sn(cst.select_leaf(cst.csa.isa[i]+1)), i) << "At position: " << i;
    }

    TEST_COUT << "Testing LCE ds" << std::endl;
    for (int i = 0; i < text.size() - 1; ++i)
    {
        EXPECT_EQ(lce_ds.lce(i, i + 1), cst.depth(cst.lca(cst.select_leaf(cst.csa.isa[i] + 1), cst.select_leaf(cst.csa.isa[i + 1] + 1)))) <<  "At position: " << i;
    }

  
}

TEST(threshold_construct_test, paper_example)
{
    std::vector<char> text = {'G', 'A', 'T', 'T', 'A', 'C', 'A', 'T', '#',
                              'G', 'A', 'T', 'A', 'C', 'A', 'T', '#',
                              'G', 'A', 'T', 'T', 'A', 'G', 'A', 'T', 'A', '#', '#'};
    std::vector<std::string> dict{"#GATTAC", "ACAT#", "AGATA##", "T#GATAC", "T#GATTAG"};
    std::vector<uint32_t> parse{1, 2, 4, 2, 5, 3, 0};
    std::vector<uint32_t> indices{0, 1, 2, 3, 4};
    std::vector<uint8_t> dict2 = {'#', '#', 'G', 'A', 'T', 'T', 'A', 'C', EndOfWord,
                                  'A', 'C', 'A', 'T', '#', EndOfWord,
                                  'A', 'G', 'A', 'T', 'A', '#', '#', EndOfWord,
                                  'T', '#', 'G', 'A', 'T', 'A', 'C', EndOfWord,
                                  'T', '#', 'G', 'A', 'T', 'T', 'A', 'G', EndOfWord, EndOfDict};
    sdsl::bit_vector b_d = {1, 0, 0, 0, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    const size_t n_phrases = 5;
    const std::vector<size_t> phrase_length{8, 5, 7, 7, 8};
    sdsl::bit_vector b_bwt = {1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1};

    const std::vector<std::pair<size_t, std::pair<size_t, size_t>>> M{
        {2, {0, 0}}, {6, {2, 2}}, {7, {3, 3}}, {7, {4, 4}}, {3, {0, 0}}, {2, {2, 3}}, {2, {4, 4}}, {3, {1, 1}}, {5, {0, 0}}, {4, {2, 2}}, {5, {3, 3}}, {5, {4, 4}}, {4, {1, 1}}, {6, {0, 0}}, {5, {2, 2}}, {6, {3, 3}}, {6, {4, 4}}, {2, {1, 1}}, {4, {0, 0}}, {3, {2, 3}}, {3, {4, 4}}, {4, {3, 3}}, {4, {4, 4}}};
    const std::vector<uint32_t> bwt_p = {3, 1, 4, 5, 2, 2};

    // Extract the reverse of the phrases
    std::vector<std::pair<std::string, uint32_t>> rev_dict;
    std::vector<char> bwt = {'A', 'T', '#', 'T', 'T', 'T', 'T', 'T', 'C', 'C',
                              'G', 'G', 'G', 'G', 'A', 'A', 'A', '#',
                              '#', '#', 'A', 'A', 'A', 'T', 'A', 'T', 'A', 'A'};
    int i = 0;
    int rank = 0;
    while (i < dict2.size() - 1)
    {
        std::string s;
        while (i < dict2.size() - 1 && dict2[i] != EndOfWord)
            s.append(1, dict2[i++]);
        i++;
        reverse(s.begin(), s.end());
        rev_dict.push_back({s, rank++});
    }
    std::sort(rev_dict.begin(), rev_dict.end());

    std::vector<uint32_t> frequencies{0, 1, 2, 1, 1, 1};
    size_t w = 2;

    TEST_COUT << "Begin paper test" << std::endl;

    TEST_COUT << "Building Pfp" << std::endl;
    pf_parsing pf(dict2, parse, frequencies, w);


    // TEST sa_ds
    std::vector<char> tmp_text = {'#', 'G', 'A', 'T', 'T', 'A', 'C', 'A', 'T', '#',
                                  'G', 'A', 'T', 'A', 'C', 'A', 'T', '#',
                                  'G', 'A', 'T', 'T', 'A', 'G', 'A', 'T', 'A', 0};

    uint8_t num_bytes = 1;
    // build cst of the Text
    pfp_lce_support lce_ds(pf);
    pfp_sa_support sa_ds(pf);

    // build cst of the Text
    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
    sdsl::construct_im(cst, static_cast<const char *>(&tmp_text[0]), num_bytes);


    TEST_COUT << "Testing LCE ds" << std::endl;
    for (int i = 0; i < tmp_text.size() - 1; ++i)
    {
        EXPECT_EQ(lce_ds.lce(i, i + 1), cst.depth(cst.lca(cst.select_leaf(cst.csa.isa[i] + 1), cst.select_leaf(cst.csa.isa[i + 1] + 1)))) << "At position: " << i;
    }
    EXPECT_EQ(pf.n, tmp_text.size());

    std::vector<uint32_t> sa(text.size(), 0);
    std::iota(sa.begin(), sa.end(), 0);
    auto cyclic_sort = [&](const size_t a, const size_t b) {
        const auto max_cmp = text.size();
        for (size_t i = 0; i < max_cmp; ++i)
        {
            if (text[(a + i) % text.size()] != text[(b + i) % text.size()])
                return text[(a + i) % text.size()] < text[(b + i) % text.size()];
        }
    };
    std::sort(sa.begin(), sa.end(), cyclic_sort);

    TEST_COUT << "Testing SA ds" << std::endl;
    for (int i = 0; i < tmp_text.size(); ++i)
    {
        EXPECT_EQ(sa_ds.sa(i), (cst.csa[i] + (cst.csa.size()) - w + 1) % (cst.csa.size())) << "At positions: " << i;
        EXPECT_EQ(text[sa_ds.sa(i)], tmp_text[cst.csa[i]]) << "At positions: " << i;
        EXPECT_EQ((sa_ds.sa(i) + pf.w - 1) % pf.n, (cst.csa[i] )% tmp_text.size()) << "At positions: " << i;
        EXPECT_EQ((sa_ds.sa(i) + pf.w - 2) % pf.n, (cst.csa[i] + tmp_text.size() - 1) % tmp_text.size()) << "At positions: " << i;
        EXPECT_EQ(sa_ds.sa(i) , sa[i]) << "At positions: " << i;
    }

    // TEST_COUT << "Testing text access" << std::endl;
    // for(size_t i = 0; i < text.size(); ++i){
    //     EXPECT_EQ((sa_ds.sa(i) + pf.w - 1) % pf.n, (sa_ds.sa(i) - 1 + pf.w) % pf.n) << "At position: " << i; // suffix number
    //     EXPECT_EQ(sa_ds.sa(i) + pf.w - 1, sa_ds.sa(i) - 1 + pf.w) << "At position: " << i; // suffix number
    // }

    TEST_COUT << "Testing text access" << std::endl;
    for (size_t i = 0; i < text.size(); ++i)
    {
        auto sn = (sa_ds.sa(i) + pf.w - 1) % pf.n; // suffix number
        auto p_i = pf.rank_b_p(sn + 1);            // phrase number
        auto id_p_i = pf.pars.p[p_i - 1];          // phrase_id of the phrase that i belongs to.
        size_t occ_in_p_i_in_D = pf.dict.select_b_d(id_p_i) + (sn - pf.select_b_p(p_i));
        auto c = pf.dict.d[occ_in_p_i_in_D];

        EXPECT_EQ(c, text[(sa_ds.sa(i) + pf.n - 1) % pf.n]) << "At position: " << i << " SN: " << sn << " p_i: " << p_i << " id_p_i: " << id_p_i << " occ_in_p_i_in_D: " << occ_in_p_i_in_D;
    }
    
    TEST_COUT << "Testing sdsl BWT" << std::endl;
    for(size_t i = 0; i < text.size(); ++i){
        EXPECT_EQ((cst.csa.bwt[i] == '\0' ? '#' : cst.csa.bwt[i]), bwt[i]) << "At position: " << i;
    }

    TEST_COUT << "Testing BWT from pf" << std::endl;
    for(size_t i = 0; i < sa_ds.size(); ++i){
        auto sn = (sa_ds.sa(i) + pf.w -1) % pf.n; // suffix number
        auto p_i = pf.rank_b_p(sn + 1);            // phrase number
        auto id_p_i = pf.pars.p[p_i - 1];          // phrase_id of the phrase that i belongs to.
        size_t occ_in_p_i_in_D = pf.dict.select_b_d(id_p_i) + (sn - pf.select_b_p(p_i));
        auto c = pf.dict.d[occ_in_p_i_in_D];

        EXPECT_EQ(c, bwt[i]) << "At position: " << i << " SN: " << sn << " p_i: " << p_i << " id_p_i: " << id_p_i << " occ_in_p_i_in_D: " << occ_in_p_i_in_D;
    }
    
    
    TEST_COUT << "Computing thresholds pfp" << std::endl;

    // This code gets timed
    std::vector<size_t> thresholds_pfp;
    {
        std::vector<uint64_t> last_seen(256, 0);
        std::vector<bool> never_seen(256, true);
        // get the charater in position i
        size_t i = 0;
        auto sn = (sa_ds.sa(i) + pf.w - 1) % pf.n;   // suffix number
        auto p_i = pf.rank_b_p(sn + 1);              // phrase number
        auto id_p_i = pf.pars.p[p_i - 1];            // phrase_id of the phrase that i belongs to.
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
                thresholds_pfp.push_back(start);
            }

            auto next_c = c;
            // Skipp the run
            while (i < sa_ds.size() && c == next_c)
            {
                ++i;
                // get the charater in position i
                auto sn = (sa_ds.sa(i) + pf.w - 1) % pf.n; // suffix number
                auto p_i = pf.rank_b_p(sn + 1);   // phrase number
                auto id_p_i = pf.pars.p[p_i - 1]; // phrase_id of the phrase that s belongs to.
                size_t occ_in_p_i_in_D = pf.dict.select_b_d(id_p_i) + (sn - pf.select_b_p(p_i));
                next_c = pf.dict.d[occ_in_p_i_in_D];
            }
            last_seen[c] = i - 1;
            never_seen[c] = false;
            c = next_c;
        }
    }

    TEST_COUT << "Computing thresholds sdsl" << std::endl;

    // This code gets timed
    std::vector<size_t> thresholds_sdsl;
    {
        std::vector<uint64_t> last_seen(256, 0);
        std::vector<bool> never_seen(256, true);

        size_t i = 0;
        // get the charater in position i
        auto c = (cst.csa.bwt[i]=='\0'?'#':cst.csa.bwt[i]);

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
                thresholds_sdsl.push_back(start);
            }

            // Skipp the run
            while (i < cst.csa.size() && c == (cst.csa.bwt[i] == '\0' ? '#' : cst.csa.bwt[i]))
                ++i;

            last_seen[c] = i - 1;
            never_seen[c] = false;
            c = (cst.csa.bwt[i] == '\0' ? '#' : cst.csa.bwt[i]);
        }
    }
    TEST_COUT << "Testing thresholds" << std::endl;
    EXPECT_EQ(thresholds_pfp.size(), thresholds_sdsl.size());
    for (int i = 0; i < std::min(thresholds_pfp.size(), thresholds_sdsl.size()); ++i)
    {
        EXPECT_EQ(thresholds_pfp[i], thresholds_sdsl[i]) << "At position: " << i;
    }
}

// TEST(lce_construct_test,lce_construct){

//     size_t w = 10;
//     pf_parsing pf(test_file, w);
//     pfp_lce_support lce_ds(pf);
//     pfp_sa_support sa_ds(pf);

//     // TEST lce_ds
//     std::vector<char> text;
//     read_fasta_file(test_file.c_str(), text);
//     std::vector<char> tmp(w-1, 2);
//     text.insert(text.begin(), tmp.begin(), tmp.end());
//     text.push_back(0);

//     TEST_COUT <<"Text size: "<< text.size() << std::endl;

//     uint8_t num_bytes = 1;
//     // build cst of the Text
//     sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
//     sdsl::construct_im(cst, static_cast<const char *>(&text[0]), num_bytes);

//     // TEST_COUT <<"Testing LCP ds" << std::endl;
//     // for (int i = 1; i < text.size() - 1; ++i)
//     // {
//     //     EXPECT_EQ(lce_ds.lce(sa_ds.sa(i), sa_ds.sa(i + 1)), cst.lcp[i + 1]);
//     // }
//     // TEST_COUT <<"Testing LCP ds" << std::endl;
//     // for (int i = 1; i < text.size() - 1; ++i)
//     // {
//     //     EXPECT_EQ(lce_ds.lce(sa_ds.sa(i), sa_ds.sa(i + 1)), cst.depth(cst.lca(cst.select_leaf(i+1), cst.select_leaf(i+2))));
//     // }

//     // TEST_COUT << "Testing SA ISA ds" << std::endl;
//     // for (int i = 0; i < text.size() - 1; ++i)
//     // {
//     //     EXPECT_EQ(cst.csa[cst.csa.isa[i]], i) << "At position: " << i;
//     //     EXPECT_EQ(cst.sn(cst.select_leaf(cst.csa.isa[i] + 1)), i) << "At position: " << i;
//     // }
//     size_t n = text.size();

//     TEST_COUT <<"Testing LCE ds" << std::endl;
//     for (int i = 0; i < text.size() - 2; ++i)
//     {
//         size_t coord_i = (i+w-1)%n;
//         size_t coord_j = (i+w)%n;

//         EXPECT_EQ(lce_ds.lce(i, i + 1), cst.depth(cst.lca(cst.select_leaf(cst.csa.isa[coord_i]+1), cst.select_leaf(cst.csa.isa[coord_j]+1)))) << "At position: " 
//                     << i << " " << coord_i << " " << i+1 << " " << coord_j;
//     }

// }













// TEST(lce_construct_test,thresholds_construct){

//     size_t w = 10;
//     TEST_COUT <<"Building pfp data structures" << std::endl;
//     pf_parsing pf(test_file, w);
//     pfp_lce_support pf_lce(pf);
//     pfp_sa_support pf_sa(pf);

//     // TEST lce_ds
//     std::vector<char> text;
//     read_fasta_file(test_file.c_str(), text);
//     std::vector<char> tmp(w-1, 2);
//     text.insert(text.begin(), tmp.begin(), tmp.end());
//     text.push_back(0);

//     TEST_COUT <<"Text size: "<< text.size() << std::endl;
//     TEST_COUT <<"Building sdsl cst" << std::endl;

//     uint8_t num_bytes = 1;
//     // build cst of the Text
//     sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
//     sdsl::construct_im(cst, static_cast<const char *>(&text[0]), num_bytes);

//     size_t n = text.size();

//     TEST_COUT <<"Computing thresholds pfp" << std::endl;

//     // This code gets timed
//     std::vector<size_t> thresholds_pfp;
//     {
//         std::vector<uint64_t> last_seen(256, 0);
//         std::vector<bool> never_seen(256, true);
//         // get the charater in position i
//         size_t i = 0;
//         auto sn = (pf_sa.sa(i) - 1 + pf.w) % pf.n;                           // suffix number
//         auto p_i = pf.rank_b_p(sn + 1); // phrase number
//         auto id_p_i = pf.pars.p[p_i - 1]; // phrase_id of the phrase that i belongs to.
//         size_t occ_in_p_i_in_D = pf.dict.select_b_d(id_p_i) + (sn - pf.select_b_p(p_i));
//         auto c = pf.dict.d[occ_in_p_i_in_D];

//         while (i < pf_sa.size())
//         {
//             // If it is not the first run
//             if (!never_seen[c])
//             {
//                 // fint the beginning of the next run of the same character.
//                 // binary search the threshold
//                 size_t start = last_seen[c], end = i, mid;
//                 while (end > start + 1)
//                 {
//                     mid = (start + end) >> 1;
//                     // LCP(start,mid)
//                     size_t lce_start_mid = pf_lce.lce(pf_sa.sa(start), pf_sa.sa(mid));
//                     // LCP(mid +1, end)
//                     size_t lce_mid_end = pf_lce.lce(pf_sa.sa(mid + 1), pf_sa.sa(end));
//                     if (lce_start_mid > lce_mid_end)
//                         start = mid;
//                     else
//                         end = mid;
//                 }
//                 thresholds_pfp.push_back(start);
//             }

//             auto next_c = c;
//             // Skipp the run
//             while (i < pf_sa.size() && c == next_c)
//             {
//                 ++i;
//                 // get the charater in position i
//                 auto sn = pf_sa.sa(i);            // suffix number
//                 auto p_i = pf.rank_b_p(sn + 1);   // phrase number
//                 auto id_p_i = pf.pars.p[p_i - 1]; // phrase_id of the phrase that s belongs to.
//                 size_t occ_in_p_i_in_D = pf.dict.select_b_d(id_p_i) + (sn - pf.select_b_p(p_i));
//                 next_c = pf.dict.d[occ_in_p_i_in_D];
//             }
//             last_seen[c] = i - 1;
//             never_seen[c] = false;
//             c = next_c;
//         }
//     }

//     TEST_COUT <<"Computing thresholds sdsl" << std::endl;

//     // This code gets timed
//     std::vector<size_t> thresholds_sdsl;
//     {
//         std::vector<uint64_t> last_seen(256, 0);
//         std::vector<bool> never_seen(256, true);

//         size_t i = 0;
//         // get the charater in position i
//         auto c = cst.csa.bwt[i];

//         while (i < cst.csa.size())
//         {



//             // If it is not the first run
//             if (!never_seen[c])
//             {
//                 // fint the beginning of the previous run of the same character.
//                 // binary search the threshold
//                 size_t start = last_seen[c], end = i, mid;
//                 while (end > start + 1)
//                 {
//                     mid = (start + end) >> 1;
//                     // LCP(start,mid)
//                     size_t lce_start_mid = cst.depth(cst.lca(cst.select_leaf(start + 1), cst.select_leaf(mid + 1)));
//                     // LCP(mid +1, end)
//                     size_t lce_mid_end = cst.depth(cst.lca(cst.select_leaf(mid + 1 + 1), cst.select_leaf(end + 1)));
//                     if (lce_start_mid > lce_mid_end)
//                         start = mid;
//                     else
//                         end = mid;
//                 }
//                 thresholds_sdsl.push_back(start);
//             }

//             // Skipp the run
//             while (i < cst.csa.size() && c == cst.csa.bwt[i])
//                 ++i;

//             last_seen[c] = i - 1;
//             never_seen[c] = false;
//             c = cst.csa.bwt[i];
//         }
//     }
//     TEST_COUT << "Testing thresholds" << std::endl;
//     EXPECT_EQ(thresholds_pfp.size(),thresholds_sdsl.size());
//     for (int i = 0; i < std::min(thresholds_pfp.size(), thresholds_sdsl.size()); ++i){
//         EXPECT_EQ(thresholds_pfp[i],thresholds_sdsl[i])<< "At position: " << i;
//     }
// }



int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " test_file " << std::endl;
        std::cout << " (1) Generates the SA, ISA and LCP;" << std::endl;
        std::cout << " (2) Generates LCE data structure and checks the result." << std::endl;
        return 1;
    }
    test_file = argv[1];
   
    return RUN_ALL_TESTS();
}


