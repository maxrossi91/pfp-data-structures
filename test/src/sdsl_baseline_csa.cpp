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

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>

#include <malloc_count.h>

int main(int argc, char const *argv[]) {

    if(argc < 2)
        error("input file required");

    std::string filename = argv[1];

    // TEST lca_ds
    verbose("Reading text from file");
    std::vector<char> text;
    read_fasta_file(filename.c_str(), text);
    verbose("Text size: " , text.size());

    uint8_t num_bytes = 1;
    // build cst of the Text
    verbose("Computing CSA of the text");
    sdsl::cache_config cc(false); // do not delete temp files after csa construction
    sdsl::csa_wt<> fm_index;
    elapsed_time(
     sdsl::construct_im(fm_index, static_cast<const char*>(&text[0]), num_bytes);
    );

    verbose("Computing LCP of the text");
    cc.delete_files = true; // delete temp files after lcp construction
    sdsl::lcp_wt<> lcp;
    elapsed_time(
     sdsl::construct_im(lcp, static_cast<const char*>(&text[0]), num_bytes);
    );

    verbose("Memory peak: ", malloc_count_peak());

  return 0;

}
