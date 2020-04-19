/* pfp - prefix free parsing sa data structure test
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
   \file sa_test.cpp
   \brief sa_test.cpp build and test prefix-free parsing sa data structure.
   \author Massimiliano Rossi
   \date 03/04/2020
*/


#include<iostream>
#include<vector>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>

#include <common.hpp>
#include <gtest/gtest.h>
#include <pfp.hpp>
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

TEST(sa_construct_test, paper_example)
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
    sdsl::bit_vector b_bwt = {1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,1,1,0,1,1,1};

    const std::vector<std::pair<size_t, std::pair<size_t, size_t>>> M{
        {2, {0, 0}}, {6, {2, 2}}, {7, {3, 3}}, {7, {4, 4}}, {3, {0, 0}},
        {2, {2, 3}}, {2, {4, 4}}, {3, {1, 1}}, {5, {0, 0}}, {4, {2, 2}},
        {5, {3, 3}}, {5, {4, 4}}, {4, {1, 1}}, {6, {0, 0}}, {5, {2, 2}},
        {6, {3, 3}}, {6, {4, 4}}, {2, {1, 1}}, {4, {0, 0}}, {3, {2, 3}},
        {3, {4, 4}}, {4, {3, 3}}, {4, {4, 4}}
    };
    const std::vector<uint32_t> bwt_p = {3, 1, 4, 5, 2, 2};

    // Extract the reverse of the phrases
    std::vector<std::pair<std::string,uint32_t>> rev_dict;
    int i = 0;
    int rank = 0;
    while(i < dict2.size()-1){
        std::string s;
        while(i < dict2.size()-1 && dict2[i] != EndOfWord) s.append(1,dict2[i++]);
        i++;
        reverse(s.begin(), s.end());
        rev_dict.push_back({s,rank++});
    }
    std::sort(rev_dict.begin(),rev_dict.end());
    
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
    
    // TEST b_bwt
    for (size_t i = 0; i < pf.b_bwt.size(); ++i)
    {
        EXPECT_EQ(pf.b_bwt[i], b_bwt[i]) << "at position: " << i;
    }
    TEST_COUT << "Test b_bwt" << std::endl;
    
    // TEST M
    for (size_t i = 0; i < pf.M.size(); ++i)
    {
        EXPECT_EQ(pf.M[i].len, M[i].first) << "at position: " << i;
        EXPECT_EQ(pf.M[i].left, M[i].second.first) << "at position: " << i;
        EXPECT_EQ(pf.M[i].right, M[i].second.second) << "at position: " << i;
    }
    TEST_COUT << "Test M" << std::endl;

    // TEST colex_document_array_helper
    std::vector<uint_t> colex_id(pf.dict.n_phrases());
    std::vector<uint_t> inv_colex_id(pf.dict.n_phrases()); // I am using it as starting positions
    for (int i = 0, j = 0; i < pf.dict.d.size(); ++i)
        if (pf.dict.d[i + 1] == EndOfWord)
        {
            colex_id[j] = j;
            inv_colex_id[j++] = i;
        }

    pf.dict.colex_document_array_helper(inv_colex_id, colex_id, 0, pf.dict.n_phrases());
    for (size_t i = 0; i < colex_id.size(); ++i){
        EXPECT_EQ(colex_id[i], rev_dict[i].second) << "at position: " << i;
    }
    TEST_COUT << "Test colex_document_array_helper" << std::endl;

    std::vector<uint32_t> sa(text.size(), 0);
    std::iota(sa.begin(), sa.end(), 0);
    auto cyclic_sort = [&](const size_t a, const size_t b) {
        const auto max_cmp = text.size();
        for (size_t i = 0; i < max_cmp; ++i) {
            if (text[(a + i) % text.size()] != text[(b + i) % text.size()])
                return text[(a + i) % text.size()] < text[(b + i) % text.size()];
        }
    };
    std::sort(sa.begin(), sa.end(), cyclic_sort);

    // TEST BWT(P) - wavelet tree
    for (size_t i = 0; i < pf.w_wt.size(); ++i)
    {
        EXPECT_EQ(pf.w_wt[i], bwt_p[i]) << "at position: " << i;
    }
    TEST_COUT << "Test BWT(P)" << std::endl;

    TEST_COUT << "Storing the PFP to file" << std::endl;
    std::string outfile = "paper_example.pf.ds";
    sdsl::store_to_file(pf, outfile.c_str());



    TEST_COUT << "Begin paper test for p1" << std::endl;

    pf_parsing pf1;
    sdsl::load_from_file(pf1, outfile);
    TEST_COUT << "Pfp built" << std::endl;

    // TEST n
    EXPECT_EQ(pf1.n, text.size());
    TEST_COUT << "Test n" << std::endl;

    // TEST n_phrases
    EXPECT_EQ(pf1.dict.n_phrases(), n_phrases);
    TEST_COUT << "Test n_phrases" << std::endl;

    // TEST b_d
    for (size_t i = 0; i < pf1.dict.b_d.size(); ++i)
    {
        EXPECT_EQ(pf1.dict.b_d[i], b_d[i]) << "at position: " << i;
    }
    TEST_COUT << "Test b_d" << std::endl;

    // TEST phrase_length
    for (size_t i = 0; i < phrase_length.size(); ++i)
    {
        EXPECT_EQ(pf1.dict.length_of_phrase(i + 1), phrase_length[i]);
    }
    TEST_COUT << "Test phrase_length" << std::endl;

    // TEST b_bwt
    for (size_t i = 0; i < pf1.b_bwt.size(); ++i)
    {
        EXPECT_EQ(pf1.b_bwt[i], b_bwt[i]) << "at position: " << i;
    }
    TEST_COUT << "Test b_bwt" << std::endl;

    // TEST M
    for (size_t i = 0; i < pf1.M.size(); ++i)
    {
        EXPECT_EQ(pf1.M[i].len, M[i].first) << "at position: " << i;
        EXPECT_EQ(pf1.M[i].left, M[i].second.first) << "at position: " << i;
        EXPECT_EQ(pf1.M[i].right, M[i].second.second) << "at position: " << i;
    }
    TEST_COUT << "Test M" << std::endl;

    // TEST colex_document_array_helper
    std::vector<uint_t> colex_id1(pf1.dict.n_phrases());
    std::vector<uint_t> inv_colex_id1(pf1.dict.n_phrases()); // I am using it as starting positions
    for (int i = 0, j = 0; i < pf1.dict.d.size(); ++i)
        if (pf.dict.d[i + 1] == EndOfWord)
        {
            colex_id1[j] = j;
            inv_colex_id1[j++] = i;
        }

    pf1.dict.colex_document_array_helper(inv_colex_id1, colex_id1, 0, pf1.dict.n_phrases());
    for (size_t i = 0; i < colex_id1.size(); ++i)
    {
        EXPECT_EQ(colex_id1[i], rev_dict[i].second) << "at position: " << i;
    }
    TEST_COUT << "Test colex_document_array_helper" << std::endl;

    // TEST BWT(P) - wavelet tree
    for (size_t i = 0; i < pf1.w_wt.size(); ++i)
    {
        EXPECT_EQ(pf1.w_wt[i], bwt_p[i]) << "at position: " << i;
    }
    TEST_COUT << "Test BWT(P)" << std::endl;



    pfp_lce_support lce_pf(pf);
    pfp_lce_support lce_pf1(pf1);

    for (int i = 0; i < text.size() - 1; ++i)
    {
        EXPECT_EQ(lce_pf.lce(i, i + 1), lce_pf1.lce(i, i + 1));
    }
    TEST_COUT << "Test LCE ds" << std::endl;

    pfp_sa_support sa_pf(pf);
    pfp_sa_support sa_pf1(pf1);

    for (int i = 0; i < sa_pf.size(); ++i)
    {
        EXPECT_EQ(sa_pf.sa(i), sa_pf1.sa(i));
    }
    TEST_COUT << "Test SA ds" << std::endl;
}

TEST(sa_construct_test,sa_construct){

    size_t w = 10;
    pf_parsing pf(test_file, w);
    pfp_sa_support sa_ds(pf);

    // TEST sa_ds
    std::vector<char> text;
    read_fasta_file(test_file.c_str(), text);
    std::vector<char> tmp(w-1,'#');
    text.insert(text.begin(),tmp.begin(),tmp.end());
    text.push_back(0);

    uint8_t num_bytes = 1;
    // build cst of the Text
    TEST_COUT << "Computing CSA of the text" << std::endl;
    sdsl::csa_wt<> csa;
    sdsl::construct_im(csa, static_cast<const char *>(&text[0]), num_bytes);

    TEST_COUT << "Testing SA ds" << std::endl;
    for (int i = 0; i < text.size(); ++i)
    {
        EXPECT_EQ(sa_ds.sa(i), (csa[i] + (text.size()) - w + 1) % (text.size())) << "At positions: " << i;
    }
}

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


