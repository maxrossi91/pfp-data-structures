/* pfp - prefix free parsing lce support
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
   \file lce_support.hpp
   \brief lce_support.hpp define and build the prefix-free parsing lce support data structures.
   \author Massimiliano Rossi
   \date 03/04/2020
*/


#ifndef _PFP_LCE_SUPPORT_HH
#define _PFP_LCE_SUPPORT_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include<pfp.hpp>

class pfp_lce_support{
protected:
  pf_parsing& pfp;

  sdsl::bit_vector b_p; // Starting sdsl::bit_vector::rank_1_type rank_b_p;rank_b_p::bit_vector::select_1_type select_b_p;
  sdsl::bit_vector::rank_1_type rank_b_p;
  sdsl::bit_vector::select_1_type select_b_p;

public:


  // This has to be changed using pfp_dictionary and pfp_parse
  pfp_lce_support(pf_parsing& pfp_):
                    pfp(pfp_),
                    b_p(pfp.n,0)
  { //

    // Build the bitvector storing the position of the beginning of each phrase.
    b_p.resize(pfp.n); // all should be initialized at false by sdsl
    for(size_t i = 0; i < 3*64; ++i) b_p[i] = false; // bug in resize
    b_p[0] = true; // phrase_0 becomes phrase 1
    size_t i = 0;
    for(int j = 0; j < pfp.pars.p.size()-1; ++j){
      // p[i]: phrase_id
      assert(pfp.pars.p[j] != 0);
      // phrase_length: select_b_d(p[i]+1)-select_b_d(p[i]);
      i += pfp.dict.length_of_phrase(pfp.pars.p[j]) - pfp.w;
      b_p[i] = true;
    }
//    b_p[n] = true; // We set the last phrase beginning in the last w characters of the string


    // Build rank and select on Sp
    rank_b_p = sdsl::bit_vector::rank_1_type(&b_p);
    select_b_p = sdsl::bit_vector::select_1_type(&b_p);
  }

  // return the longest common prefix of suffix i and j of T
  size_t lce(size_t i, size_t j){
    if( i == j )
      return pfp.n-i;

    // find the phrases of which i and j belongs to
    size_t p_i = rank_b_p(i+1);// - 1; // rank of the phrase in T that i belongs to. Span: [1..|P|]
    size_t p_j = rank_b_p(j+1);// - 1; // rank of the phrase in T that j belongs to. Span: [1..|P|]
    // NOTE: i+1 since rank return the number of 1s in v[0..i-1]


    // NOTE: the -1 is because the bitvector p start with a 1 but the first phrase is in position 0 in p
    auto id_p_i = pfp.pars.p[p_i - 1]; // phrase_id of the phrase that i belongs to.
    auto id_p_j = pfp.pars.p[p_j - 1]; // phrase_id of the phrase that j belongs to.

    auto tmp_s_pi = select_b_p(p_i);
    auto tmp_s_pj = select_b_p(p_j);
    // find the length of the suffixes iof id_p_i and id_p_j starting at i and j
    // Length of the phrrase - length of the prefix.
    size_t len_suff_i_in_p_i = pfp.dict.length_of_phrase(id_p_i) - (i - select_b_p(p_i) + 1);
    size_t len_suff_j_in_p_j = pfp.dict.length_of_phrase(id_p_j) - (j - select_b_p(p_j) + 1);
    size_t k = std::min(len_suff_i_in_p_i,len_suff_j_in_p_j);

     // find the occurrence of i and j in the concatenation of the phrases of D
     // Starting position of the phrase in D + length of the prefix.
     size_t occ_in_p_i_in_D = pfp.dict.select_b_d(id_p_i) + (i - select_b_p(p_i) + 1);
     size_t occ_in_p_j_in_D = pfp.dict.select_b_d(id_p_j) + (j - select_b_p(p_j) + 1);

     // compute the LCP of the suffix up to the phrase boundaries.
     auto lcp_ppi_ppj = 0;
     if(occ_in_p_i_in_D == occ_in_p_j_in_D){
       // if two positions map to the each position in the same phrase,
       // we have to ignore the lcp in D because the lcp is the length of the suffix of the phrase.
       lcp_ppi_ppj = k;
     }else{
       auto lcp_left = std::min(pfp.dict.isaD[occ_in_p_i_in_D],pfp.dict.isaD[occ_in_p_j_in_D])+1;
       auto lcp_right = max(pfp.dict.isaD[occ_in_p_i_in_D],pfp.dict.isaD[occ_in_p_j_in_D]);
       size_t lcp_ppi_ppj_i = pfp.dict.rmq_lcp_D(lcp_left,lcp_right);
       lcp_ppi_ppj = pfp.dict.lcpD[lcp_ppi_ppj_i];
     }

     if(lcp_ppi_ppj < k){ // if the lcp of the suffixes up to the phrase boundaries is shorter than the shortest suffix length
       return lcp_ppi_ppj;
     }else{ // otherwise, find the first phrase where the suffixes differ

       k -= pfp.w; // this align the lcp values
       // Check if one of those is the last should be dobe by the $ in P
       auto p_ip1_in_sa = pfp.pars.isaP[p_i]; // isaP is 0-based
       auto p_jp1_in_sa = pfp.pars.isaP[p_j];

       auto lcp_left = std::min(p_ip1_in_sa,p_jp1_in_sa)+1;
       auto lcp_right = max(p_ip1_in_sa,p_jp1_in_sa);
       size_t lcp_pi_pj_i = pfp.pars.rmq_lcp_P(lcp_left,lcp_right);
       auto lcp_pi_pj = pfp.pars.lcpP[lcp_pi_pj_i];

       // size_t l_com_phrases = 0
       // for l in range(lcp_pi_pj):
       // #               print(p_i + 1 + l)
       // l_com_phrases += len(self.D[self.P[p_i + 1 + l]]) - w

       size_t l_com_phrases = select_b_p(p_i + 1 + lcp_pi_pj) - select_b_p(p_i + 1); // Check if there is some error it might be here


       auto a = pfp.pars.p[p_i + lcp_pi_pj]; // p is 0-based
       auto b = pfp.pars.p[p_j + lcp_pi_pj];

       assert (a != 0 and b != 0);
       // if a == "#" or b == "#": # This cannot be possible due to the last $ in T'
       // return k + l_com_phrases
       // else:
       // Compute the lcp between phrases a and b
       auto a_in_sa = pfp.dict.isaD[pfp.dict.select_b_d(a)]; // position of the phrase a in saD
       auto b_in_sa = pfp.dict.isaD[pfp.dict.select_b_d(b)]; // position of the phrase b in saD

       lcp_left = std::min(a_in_sa,b_in_sa) + 1;
       lcp_right = max(a_in_sa,b_in_sa);

       size_t lcp_a_b_i = pfp.dict.rmq_lcp_D(lcp_left, lcp_right);
       auto lcp_a_b = pfp.dict.lcpD[lcp_a_b_i];

       return k + l_com_phrases + lcp_a_b;
     }



  }

};

#endif /* end of include guard: _PFP_LCE_SUPPORT_HH */
