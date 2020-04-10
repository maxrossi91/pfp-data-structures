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

#include <pfp.hpp>
#include <wt.hpp>
#include <lce_support.hpp>

class pfp_sa_support {
protected:
  pf_parsing& pfp;
  pfp_lce_support& lce_support;
public:
  pfp_sa_support(pf_parsing & pfp_, pfp_lce_support & lce_support)
    : pfp(pfp_), lce_support(lce_support)
  { }

  // TODO: get type return type
  uint32_t sa(uint32_t i) {
    // b_bwt.rank(i) - 1
    // i - b_bwt.select(b_bwt.rank(i))
    // -> to find lex rank of the proper phrase suffix a of length at least w
    //    that starts at SA[i], and lex rank j of S[SA[i] .. n - 1]
    //    among suffixes of S starting with a
    pfp.w_wt.print_leafs();

    const auto rank_i = pfp.b_bwt_rank_1(i + 1);
    const auto lex_rank_i = rank_i - 1;
    // lex rank - 0-based
    const auto interval_rank = i - pfp.b_bwt_select_1(rank_i);

    std::cout << "SA[" << i << "]:\n"
              << "> rank_i: " << rank_i << std::endl
              << "> lex_rank_i: " << lex_rank_i << std::endl
              << "> interval_rank: " << interval_rank << std::endl;

    const auto & m = pfp.M[lex_rank_i];
    std::cout << "M: len = " << m.len << " | [" << m.left
              << ", " << m.right << "]" << std::endl;

    // WT have 0 delimiter, lex smallest, so it shifted intervals
    // rank + 1 because its 0-based and select is 1-based
    const auto k = pfp.w_wt.range_select(m.left, m.right, interval_rank + 1);
    uint32_t p_i = pfp.pars.saP[k + 1];
    std::cout << "k: " << k << " | p_i: " << p_i << std::endl;

    // const auto occ_k = lce_support.select_b_p(p_i); // start of phrase in S
    const auto occ_k_next = lce_support.select_b_p(p_i + 1); // start of next phrase in S

    // because b_p starts with trigger string (cyclic S):
    // - start of next phrase is select(i + 2) = next
    // - SA[i] = next - M.len - w -> trigger string at start
    // because b_p starts with trigger string - we need to -w positions

    std::cout << "occ k + 1: " << occ_k_next << std::endl;

    std::cout << "SA[" << i << "] = " << occ_k_next - (m.len - pfp.w) - pfp.w << std::endl;
    // const auto phrase_id = pfp.pars.p[p_i];
    // const auto loc = pfp.select_b_p(p_i);
  }
};


#endif /* end of include guard: _PFP_SA_SUPPORT_HH */
