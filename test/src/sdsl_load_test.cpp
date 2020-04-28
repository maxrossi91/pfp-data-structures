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
#include <sdsl/io.hpp>
#include <sdsl/suffix_trees.hpp>

#include <pfp.hpp>
#include <lce_support.hpp>
#include <sa_support.hpp>

#include <malloc_count.h>


int main(int argc, char* const argv[]) {


  Args args;
  parseArgs(argc, argv, args);


  // verbose("Window size set to: " , args.w);

  // pf_parsing pf1(args.filename, args.w);

  // verbose("Storing the PFP to file");
  // sdsl::store_to_file(pf1, outfile.c_str());

  verbose("Loading PFP data structures from:", args.filename);
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
  std::string filename = args.filename + ".sdsl.cst";
  sdsl::load_from_file(cst, filename);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("PFP DS construction complete");
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
  auto time = std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count();


  for (int i = 0; i < cst.csa.size(); ++i)
  {
    info(cst.csa[i]);
  }

  return 0;

}
