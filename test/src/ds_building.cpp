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

//********** argument options ********************
// struct containing command line parameters and other globals
struct Args{
  std::string filename = "";
  size_t w = 10;         // sliding window size and its default
  bool store = false;    // store the data structure in the file
  bool memo = false;     // print the memory usage
  // bool is_fasta = false; // read a fasta file
};

void parseArgs(int argc, char* const argv[], Args &arg){
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-w wsize] [-s store] [-m memo]\n\n" +
                    "Computes the pfp data structures of infile, provided that infile.parse, infile.dict, and infile.occ exists.\n" +
                    " store: [boolean] - store the data structure in infile.pfp.ds. (def. false)\n" +
                    "  memo: [boolean] - print the data structure memory usage. (def. false)\n" +
                    " wsize: [integer] - sliding window size (def. 10)\n");

  std::string sarg;
  while ((c = getopt(argc, argv, "w:mh")) != -1){
    switch (c){
    case 's':
      arg.store = true;
      break;
    case 'w':
      sarg.assign(optarg);
      arg.w = stoi(sarg);
      break;
    case 'm':
      arg.memo = true;
      break;
    // case 'f':
    //   arg.is_fasta = true;
    //   break;
    case 'h':
      error(usage);
    case '?':
      error("Unknown option.\n" , usage);
      exit(1);
    }
  }
  // the only input parameter is the file name
  if (argc == optind + 1){
    arg.filename.assign(argv[optind]);
  }
  else{
    error("Invalid number of arguments\n" , usage);
  }
}

//********** end argument options ********************

int main(int argc, char* const argv[]) {


  Args args;
  parseArgs(argc, argv, args);


  verbose("Window size set to: " , args.w);

  verbose("Computing PFP data structures");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  pf_parsing pf(args.filename, args.w);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("PFP DS construction complete");
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());


  verbose("Providing LCE support");
  _elapsed_time(
    pfp_lce_support lce_ds(pf)
  );

  verbose("Providing SA support");
  _elapsed_time(
    pfp_sa_support pfp_sa(pf)
  );

  verbose("Memory peak: ", malloc_count_peak());

  if(args.memo){
    verbose("Dictionary size (bytes)  : ", sdsl::size_in_bytes(pf.dict));
    verbose("Parse size (bytes)       : ", sdsl::size_in_bytes(pf.pars));
    verbose("Wavelet tree size (bytes): ", sdsl::size_in_bytes(pf.w_wt));

    verbose("PFP DS size (bytes): ", sdsl::size_in_bytes(pf));
  }

  if(args.store){
    verbose("Storing the PFP to file");
    std::string outfile = args.filename + ".pf.ds";
    sdsl::store_to_file(pf, outfile.c_str());
  }

  return 0;

}
