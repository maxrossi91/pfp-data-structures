/* sdsl-baseline_cst - prefix free parsing data structures
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
   \file sdsl_baseline_cst.cpp
   \brief sdsl_baseline_cst.cpp build an sdsl CST.
   \author Massimiliano Rossi
   \date 02/04/2020
*/

#include<iostream>

#define VERBOSE

#include <common.hpp>
#include <strdup.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/cst_sct3.hpp>

#include <malloc_count.h>

int main(int argc, char *const argv[])
{

    Args args;
    parseArgs(argc, argv, args);

    // TEST lca_ds
    verbose("Reading text from file");
    std::vector<char> text;
    read_fasta_file(args.filename.c_str(), text);
    verbose("Text size: " , text.size());

    uint8_t num_bytes = 1;
    // build cst of the Text
    verbose("Computing CST of the text");
    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
    auto time = _elapsed_time(
     sdsl::construct_im(cst, static_cast<const char*>(&text[0]), num_bytes);
    );

    auto mem_peak = malloc_count_peak();
    verbose("Memory peak: ", malloc_count_peak());

    size_t space = 0;
    if(args.memo){
        space = size_in_bytes(cst);
        verbose("CST size: ", space);
    }

    if(args.store){
        verbose("Storing the CST to file");
        std::string outfile = args.filename + ".sdsl.cst";
        store_to_file(cst, outfile.c_str());
    }
    
    if(args.csv)
        std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
}
