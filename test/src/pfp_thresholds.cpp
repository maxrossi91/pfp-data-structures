/* pfp_thresholds - Test of construction of MEM thresholds from prefix free parsing data structures
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
   \file ptp_thresholds.cpp
   \brief ptp_thresholds.cpp test construction of MEM-thresholds from prefix-free parsing data structures construction.
   \author Massimiliano Rossi
   \date 25/04/2020
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

  verbose("Loading PFP data structures from file");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  pf_parsing<> pf;
  std::string filename = args.filename + pf.filesuffix();
  sdsl::load_from_file(pf, filename);

  pfp_sa_support<> sa_ds(pf);
  pfp_lce_support<> lce_ds(pf);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("PFP DS loading complete");
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Loading BWT from file");
  std::vector<uint8_t> bwt;
  filename = args.filename + std::string(".dict");
  // Load BWT from file from file
  _elapsed_time(
    read_file(filename.c_str(), bwt);
  );

  verbose("Building the thresholds");
  
  std::vector<size_t> thresholds;
  std::vector<uint64_t> last_seen(256, 0);
  std::vector<bool> never_seen(256, true);

  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

  // This code gets timed
  // get the charater in position i
  size_t i = 0;
  auto c = bwt[i];

  while (i < sa_ds.size())
  {
    // If it is not the first run
    if (!never_seen[c])
    {
      // fint the beginning of the next run of the same character.
      // binary search the threshold
      size_t start = last_seen[c], end = i, mid;
      while (end > start + 1)
      {
        mid = (start + end) >> 1;
        // LCP(start,mid)
        size_t lce_start_mid = lce_ds.lce(sa_ds.sa(start), sa_ds.sa(mid));
        // LCP(mid +1, end)
        size_t lce_mid_end = lce_ds.lce(sa_ds.sa(mid + 1), sa_ds.sa(end));
        if (lce_start_mid > lce_mid_end)
          start = mid;
        else
          end = mid;
      }
      thresholds.push_back(start);
    }

    auto next_c = c;
    // Skipp the run
    while (i < sa_ds.size() && c == next_c)
    {
      ++i;
      // get the charater in position i
      next_c = bwt[i];
    }
    last_seen[c] = i - 1;
    never_seen[c] = false;
    c = next_c;
  }


  std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration<double, std::ratio<1>>(t_end - t_start).count();
  verbose("Elapsed time (s): ", time);

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  size_t space = 0;
  if(args.memo){
    space = thresholds.size() * sizeof(thresholds[0]);
    verbose("Thresholds size (bytes): ", space);
  }

  if(args.store){
    verbose("Storing the Thresholds to file");
    std::string outfile = args.filename + std::string(".thr");
    sdsl::store_to_file(thresholds, outfile.c_str());
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;

}
