/* pfp - prefix free parsing suffix array support
    Copyright (C) 2020 Ondřej Cvacho

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
   \file sa_support.hpp
   \brief sa_support.hpp define and build the prefix-free parsing suffix array support data structures.
   \author Ondřej Cvacho
   \date 03/04/2020
*/


#ifndef _PFP_SA_SUPPORT_HH
#define _PFP_SA_SUPPORT_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include<pfp.hpp>
#include<wt.hpp>

class pfp_sa_support {
protected:
  pf_parsing& pfp;

  pfp_wt wt;
public:
  pfp_sa_support(pf_parsing & pfp_)
    : pfp(pfp_)
  {
    verbose("Creating PFP SA data structure");
    // create alphabet (phrases)
    std::vector<uint32_t> alphabet(pfp.dict.n_phrases());
    std::iota(alphabet.begin(), alphabet.end(), 1);

    // sort alphabet based on co-lexicographical order of phrases
    // pfp.pars.p - seq of phrase ids
    // pfp.dict.b_d
    // co-lexi compare function
    auto co_lexi_dict_cmp = [&](const uint32_t i, const uint32_t j) {
      auto i_start = pfp.dict.select_b_d(i);
      auto i_end = i_start + pfp.dict.length_of_phrase(i) - 1;
      auto j_start = pfp.dict.select_b_d(j);
      auto j_end = j_start + pfp.dict.length_of_phrase(j) - 1;

      auto i_r_begin = pfp.dict.d.rend() - i_end - 1;
      auto i_r_end = pfp.dict.d.rend() - i_start;
      auto j_r_begin = pfp.dict.d.rend() - j_end - 1;
      auto j_r_end = pfp.dict.d.rend() - j_start;

      return std::lexicographical_compare(i_r_begin, i_r_end, j_r_begin, j_r_end);
    };
    std::sort(alphabet.begin(), alphabet.end(), co_lexi_dict_cmp);

    wt.construct(alphabet, pfp.pars.p);
    // wt.print_leafs();
  }
};


#endif /* end of include guard: _PFP_SA_SUPPORT_HH */
