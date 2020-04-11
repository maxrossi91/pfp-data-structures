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

  size_t size() const {
    return pfp.n;
  }

  uint32_t sa(uint32_t i) {
    // i + 1 -> rank is for (0 ... i - 1) interval
    const auto rank_i = pfp.b_bwt_rank_1(i + 1);
    const auto lex_rank_i = rank_i - 1;
    const auto interval_rank = i - pfp.b_bwt_select_1(rank_i); // lex rank - 0-based
    const auto & m = pfp.M[lex_rank_i];

    // WT have 0 delimiter, lex smallest, so it shifted intervals
    // rank + 1 because its 0-based and select is 1-based
    const auto k = pfp.w_wt.range_select(m.left, m.right, interval_rank + 1);

    uint32_t p_i;
    if (pfp.pars.saP[k + 1] > 0)
      p_i = pfp.pars.saP[k + 1] - 1;
    else
      p_i = pfp.pars.p.size() - 2;

    size_t occ_k_next;
    if (p_i + 2 > pfp.pars.p.size() - 1)
      occ_k_next = pfp.n;
    else
      occ_k_next = lce_support.select_b_p(p_i + 2); // start of next phrase in S

    if (occ_k_next < m.len)
      return pfp.n - (m.len - occ_k_next);
    else
      return occ_k_next - m.len; // b_p starts with trigger string (occ_k_next - (m.len - pfp.w) - pfp.w)
  }
};

#endif /* end of include guard: _PFP_SA_SUPPORT_HH */
