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

#include <common.hpp>
#include <gtest/gtest.h>
#include <pfp.hpp>
#include <lce_support.hpp>

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


    pfp_lce_support lce_ds(pf);

    // TEST b_p
    for (size_t i = 0; i < lce_ds.b_p.size(); ++i)
    {
        EXPECT_EQ(lce_ds.b_p[i], b_p[i]) << "at position: " << i;
    }
    TEST_COUT << "Test b_p" << std::endl;
    
    
    // TEST lce_ds
    // build lcp of the Text
    uint8_t num_bytes = 1;
    // build cst of the Text
    verbose("Computing CSA of the text");
    sdsl::cache_config cc(false); // do not delete temp files after csa construction
    sdsl::csa_wt<> csa;
    sdsl::construct_im(csa, static_cast<const char *>(&text[0]), num_bytes);

    verbose("Computing LCP of the text");
    cc.delete_files = true; // delete temp files after lcp construction
    sdsl::lcp_wt<> lcp;
    sdsl::construct_im(lcp, static_cast<const char *>(&text[0]), num_bytes);

    verbose("Testing LCE ds");
    for (int i = 1; i < text.size() - 1; ++i)
    {
        EXPECT_EQ(lce_ds.lce(csa[i], csa[i + 1]), lcp[i + 1]) << "At positions: " << csa[i] << " " << csa[i + 1];
    }
}

TEST(lce_construct_test,lce_construct){

    size_t w = 10;
    pf_parsing pf(test_file, w);
    pfp_lce_support lce_ds(pf);

    // TEST lce_ds
    std::vector<char> text;
    read_fasta_file(test_file.c_str(), text);
    verbose("Text size: ", text.size());

    uint8_t num_bytes = 1;
    // build cst of the Text
    verbose("Computing CSA of the text");
    sdsl::cache_config cc(false); // do not delete temp files after csa construction
    sdsl::csa_wt<> csa;
    sdsl::construct_im(csa, static_cast<const char *>(&text[0]), num_bytes);

    verbose("Computing LCP of the text");
    cc.delete_files = true; // delete temp files after lcp construction
    sdsl::lcp_wt<> lcp;
    sdsl::construct_im(lcp, static_cast<const char *>(&text[0]), num_bytes);

    verbose("Testing LCE ds");
    for (int i = 1; i < text.size() - 1; ++i)
    {
        auto a = lce_ds.lce(csa[i], csa[i + 1]);
        auto b = lcp[i + 1];
        EXPECT_EQ(a, b);
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


