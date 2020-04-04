/* ds_building - Test of construction of prefix free parsing data structures
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
   \file ds_building.cpp
   \brief ds_building.cpp test prefix-free parsing data structures construction.
   \author Massimiliano Rossi
   \date 20/03/2020
*/

#include<iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C" {
  #include <gsacak.h>
}
#include <malloc_count.h>

#include <pfp.hpp>
#include <lce_support.hpp>
#include <sa_support.hpp>

// TODO:
//  - Rozdeleni alphabet je potreba kontrolovat na pritomnost v te konkretni casti alphabety
void create_W_simple() {
  std::vector<std::string> dict {"##GATTAC", "ACAT#", "AGATA##", "T#GATAC", "T#GATTAG"};
  std::vector<uint32_t> parse {0, 1, 3, 1, 4, 2};
  std::vector<uint32_t> indices {0, 1, 2, 3, 4};

  std::vector<uint8_t> dict2 = {'#', '#', 'G', 'A', 'T', 'T', 'A', 'C', EndOfWord,
                                'A', 'C', 'A', 'T', '#', EndOfWord,
                                'A', 'G', 'A', 'T', 'A', '#', '#', EndOfWord,
                                'T', '#', 'G', 'A', 'T', 'A', 'C', EndOfWord,
                                'T', '#', 'G', 'A', 'T', 'T', 'A', 'G', EndOfWord
                                };
  sdsl::bit_vector b_d = {1, 0, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 0, 0, 0, 0,
                          1, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 0, 0, 0, 0, 0, 0, 0
                          };
  sdsl::bit_vector::rank_1_type rank_b_d(&b_d);
  sdsl::bit_vector::select_1_type select_b_d(&b_d);
  const size_t n_phrases = 5;
  const std::vector<size_t> phrase_length {8, 5, 7, 7, 8};

  // co-lexi compare function
  auto co_lexi_dict_cmp = [&](const uint32_t & i, const uint32_t & j) -> bool {
    // compare i-th phrase with j-th phrase
    // get start and end of phrases
    auto i_start = select_b_d(i + 1);
    auto i_end = i_start + phrase_length[i] - 1;
    auto j_start = select_b_d(j + 1);
    auto j_end = j_start + phrase_length[j] - 1;

    std::cout << "------------- i: " << i << " j: " << j << std::endl
              << "i_s: " << i_start << " | i_e: " << i_end << std::endl
              << "j_s: " << j_start << " | j_e: " << j_end << std::endl;

    auto i_r_begin = dict2.rend() - i_end - 1; // remove i EndOfWords
    auto i_r_end = dict2.rend() - i_start;
    auto j_r_begin = dict2.rend() - j_end - 1;
    auto j_r_end = dict2.rend() - j_start;
    std::cout << "> comparing" << std::endl
              << "  - ";
    std::copy(i_r_begin, i_r_end, std::ostream_iterator<uint8_t>(std::cout, ""));
    std::cout << std::endl
              << "  - ";
    std::copy(j_r_begin, j_r_end, std::ostream_iterator<uint8_t>(std::cout, ""));
    std::cout << std::endl;

    return std::lexicographical_compare(i_r_begin, i_r_end, j_r_begin, j_r_end);
  };
  std::sort(indices.begin(), indices.end(), co_lexi_dict_cmp);

  std::cout << "0: " << select_b_d(1) << " | 1: " << select_b_d(2) << std::endl;
  std::cout << "last: " << select_b_d(5) << " | last+1: " << select_b_d(6) << std::endl;

  for (const auto i : indices)
    std::cout << i << " | ";
  std::cout << std::endl;

  // create wavelet tree based on co-lexi sorted phrases
  // pfp_wt wt(indices, parse);
  // wt.print_leafs();
}

int main(int argc, char const *argv[]) {

  // create_W_simple();
  // return 0;

  if(argc < 2)
    error("input file required");

  std::string filename = argv[1];

  size_t w = 10;
  pf_parsing pf(filename,w);

  verbose("Dictionary size (bytes):       " , sizeof(pf.dict.d[0]) * pf.dict.d.size());
  verbose("Dictionary SA size (bytes):    " , sizeof(pf.dict.saD[0]) * pf.dict.saD.size());
  verbose("Dictionary ISA size (bytes):   " , sizeof(pf.dict.isaD[0]) * pf.dict.isaD.size());
  verbose("Dictionary LCP size (bytes):   " , sizeof(pf.dict.lcpD[0]) * pf.dict.lcpD.size());
  verbose("Dictionary DA size (bytes):    " , sizeof(pf.dict.daD[0]) * pf.dict.daD.size());
  verbose("Dictionary co-DA size (bytes): " , sizeof(pf.dict.colex_daD[0]) * pf.dict.colex_daD.size());
  verbose("Parsing size (bytes):          " , sizeof(pf.pars.p[0]) * pf.pars.p.size());
  verbose("Parsing SA size (bytes):       " , sizeof(pf.pars.saP[0]) * pf.pars.saP.size());
  verbose("Parsing ISA size (bytes):      " , sizeof(pf.pars.isaP[0]) * pf.pars.isaP.size());
  verbose("Parsing LCP size (bytes):      " , sizeof(pf.pars.lcpP[0]) * pf.pars.lcpP.size());


  size_t n = pf.n;

  verbose("Providing LCE support");
  _elapsed_time(
    pfp_lce_support lce_ds(pf)
  );

  verbose("Computing W");
  _elapsed_time(
    pfp_sa_support pfp_sa(pf)
  );

  // Building b_bwt

#if 0
// TEST lca_ds
std::vector<unsigned char> text;
read_fasta_file(filename.c_str(), text);
verbose("Text size: " , text.size());

// build lcp of the Text
text.push_back(0);
uint_t t_len = text.size();
std::vector<uint_t> saT(t_len); //uint_t *saT = new uint_t[t_len];
std::vector<int_t> lcpT(t_len); //int_t *lcpT = new int_t[t_len];

verbose("Computing SA, ISA, and LCP of the text");
// time_t  start = time(NULL);
sacak(&text[0],&saT[0],t_len);

// Computing isaT
std::vector<uint_t> isaT(t_len); //uint_t *isaT = new uint_t[t_len];
for(int i = 0; i < t_len; ++i){
  isaT[saT[i]] = i;
}

LCP_array(&text[0], isaT, saT, t_len, lcpT);
sdsl::rmq_succinct_sct<> rmq_lcp_T = sdsl::rmq_succinct_sct<>(&lcpT);

/*
// Compute the Dictionary
std::vector<std::string> di;

std::string s;
for(size_t i=0;i<saD.size();i++){
  if(d[i]==EndOfWord) {
    di.push_back(s);
    s.clear();
  }else{
    s.append(1,d[i]);
  }
}
*/

verbose("Testing LCE ds");
// for(int i = 1; i < text.size()-1 ; ++i){
//   for(int j = i+1; j < text.size(); ++j){
//     auto a = lce_ds.lce(saT[i], saT[j]);
//     auto b = lcpT[rmq_lcp_T(i+1,j)];//lcpT[i+1];
//     assert( a == b);
//   }
// }
for(int i = 1; i < text.size()-1 ; ++i){
    auto a = lce_ds.lce(saT[i], saT[i+1]);
    auto b = lcpT[i+1];//lcpT[rmq_lcp_T(i+1,j)];//lcpT[i+1];
    assert( a == b);
}

#endif



  return 0;

}
