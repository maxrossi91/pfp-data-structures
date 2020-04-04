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

#include <common.hpp>
#include <gtest/gtest.h>
#include <pfp.hpp>
#include <lce_support.hpp>

extern "C" {
    #include<gsacak.h>
}

std::string test_file;

TEST(lce_construct_test,lce_construct){

    size_t w = 10;
    pf_parsing pf(test_file, w);
    pfp_lce_support lce_ds(pf);

    // TEST lce_ds
    std::vector<unsigned char> text;
    read_fasta_file(test_file.c_str(), text);
    verbose("Text size: ", text.size());

    // build lcp of the Text
    text.push_back(0);
    uint_t t_len = text.size();
    std::vector<uint_t> saT(t_len); //uint_t *saT = new uint_t[t_len];
    std::vector<int_t> lcpT(t_len); //int_t *lcpT = new int_t[t_len];

    verbose("Computing SA, ISA, and LCP of the text");
    // time_t  start = time(NULL);
    sacak(&text[0], &saT[0], t_len);

    // Computing isaT
    std::vector<uint_t> isaT(t_len); //uint_t *isaT = new uint_t[t_len];
    for (int i = 0; i < t_len; ++i)
    {
        isaT[saT[i]] = i;
    }

    LCP_array(&text[0], isaT, saT, t_len, lcpT);

    verbose("Testing LCE ds");
    for (int i = 1; i < text.size() - 1; ++i)
    {
        auto a = lce_ds.lce(saT[i], saT[i + 1]);
        auto b = lcpT[i + 1];
        EXPECT_EQ(a,b);
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


