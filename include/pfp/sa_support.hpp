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
public:
  pfp_sa_support(pf_parsing & pfp_)
    : pfp(pfp_)
  {
    verbose("Creating PFP SA data structure");
  }

  // TODO: get type return type
  uint32_t sa(uint32_t idx) {
    
  }
};


#endif /* end of include guard: _PFP_SA_SUPPORT_HH */
