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
#include <strdup.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include <pfp.hpp>
#include <lce_support.hpp>
#include <sa_support.hpp>

#include <malloc_count.h>


int main(int argc, char* const argv[]) {


  Args args;
  parseArgs(argc, argv, args);


  verbose("Window size set to: " , args.w);

  verbose("Computing PFP data structures");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  pf_parsing<pfp_wt_sdsl_2> pf(args.filename, args.w);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("PFP DS construction complete");
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
  auto time = std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count();

  verbose("Providing LCE support");
  _elapsed_time(
      pfp_lce_support<pfp_wt_sdsl_2> lce_ds(pf)
  );

  verbose("Providing SA support");
  _elapsed_time(
      pfp_sa_support<pfp_wt_sdsl_2> pfp_sa(pf)
  );

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  size_t space = 0;
  if(args.memo){
    verbose("Dictionary size (bytes)  : ", sdsl::size_in_bytes(pf.dict));
    verbose("Parse size (bytes)       : ", sdsl::size_in_bytes(pf.pars));
    verbose("Wavelet tree size (bytes): ", sdsl::size_in_bytes(pf.w_wt));
    verbose("M size (bytes)           : ", sdsl::size_in_bytes(pf.M));

    space = sdsl::size_in_bytes(pf);
    verbose("PFP DS size (bytes): ", space);
  }

  if(args.store){
    verbose("Storing the PFP to file");
    std::string outfile = args.filename + pf.filesuffix();
    sdsl::store_to_file(pf, outfile.c_str());
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;

}
