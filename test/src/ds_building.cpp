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

#include <pfp.hpp>
#include <lce_support.hpp>
#include <sa_support.hpp>

#include <malloc_count.h>



int main(int argc, char const *argv[]) {

  std::string usage("usage: " + std::string(argv[0]) + " infile [-w wsize]\n\n" +
                    "Computes the pfp data structures of infile, provided that infile.parse, infile.dict, and infile.occ exists.\n" +
                    " wsize: [integer] - sliding window size (def. 10)\n");

  if(argc < 2 || argc == 3 || argc > 4)
    error(usage);

  // Parse argv
  std::string filename = argv[1];
  
  size_t w = 10;
  
  if(argc == 4){
    std::string opt(argv[2]);
    if(opt == "-w") w = std::stoi(argv[3]);
    else error(usage);
  }

  verbose("Window size set to: " , w);




  pf_parsing pf(filename,w);

  verbose("PFP DS construction complete");
  size_t n = pf.n;

  verbose("Providing LCE support");
  _elapsed_time(
    pfp_lce_support lce_ds(pf)
  );

  verbose("Providing SA support");
  _elapsed_time(
    pfp_sa_support pfp_sa(pf)
  );

  verbose("Memory peak: ", malloc_count_peak());

  verbose("PFP DS size (bytes): ", sdsl::size_in_bytes(pf));

  verbose("Storing the PFP to file");
  std::string outfile = filename + ".pf.ds";
  sdsl::store_to_file(pf, outfile.c_str());
  return 0;

}
