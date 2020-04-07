/* pfp - prefix free parsing 
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
   \file pfp.hpp
   \brief pfp.hpp define and build the prefix-free parsing data structures.
   \author Massimiliano Rossi
   \date 03/04/2020
*/


#ifndef _PFP_HH
#define _PFP_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C" {
    #include<gsacak.h>
}

#include<dictionary.hpp>
#include<parse.hpp>

class pf_parsing{
public:
  struct M_entry_t{
    uint_t len;
    uint_t left; // left and right are the extreemes of the range
    uint_t right;
  };


  dictionary dict;
  parse pars;
  std::vector<uint32_t> freq;
  size_t n; // Size of the text
  size_t w; // Size of the window

  sdsl::bit_vector b_bwt;
  std::vector<M_entry_t> M;

  pf_parsing(std::vector<uint8_t> &d_,
             std::vector<uint32_t> &p_,
             std::vector<uint32_t> &freq_,
             size_t w_) : 
            dict(d_, w_),
            pars(p_, dict.n_phrases() + 1),
            freq(freq_),
            w(w_)
  {

    // Uploading the frequency file
    assert(freq[0] == 0);

    // Compute the length of the string;
    compute_n();

    verbose("Computing b_bwt and M of the parsing");
    _elapsed_time(build_b_bwt_and_M()); 
  }

  pf_parsing( std::string filename, size_t w_):
              dict(filename, w_),
              pars(filename,dict.n_phrases()+1),
              freq(1,0),
              w(w_)
  {

    // Uploading the frequency file
    uint32_t *occ;
    size_t d_words;
    std::string tmp_filename = filename + std::string(".occ");
    read_file<uint32_t> (tmp_filename.c_str(), occ, d_words);
    freq.insert(freq.end(),occ, occ + d_words);


    // Compute the length of the string;
    compute_n();

    verbose("Computing b_bwt and M of the parsing");
    _elapsed_time(build_b_bwt_and_M()); 
  }

  void compute_n(){
    // Compute the length of the string;
    n = 0;
    for (int j = 0; j < pars.p.size() - 1; ++j)
    {
      // parse.p[j]: phrase_id
      assert(pars.p[j] != 0);
      n += dict.length_of_phrase(pars.p[j]) - w;
    }
    //n += w; // + w because n is the length including the last w markers
    //n += w - 1; // Changed after changind b_d in dict // -1 is for the first dollar + w because n is the length including the last w markers
  }

  void build_b_bwt_and_M()
  {
    // Build the bitvector storing the position of the beginning of each phrase.
    b_bwt.resize(n);
    for (size_t i = 0; i < b_bwt.size(); ++i)
      b_bwt[i] = false; // bug in resize

    assert(dict.d[dict.saD[0]] == EndOfDict);
    size_t i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
    size_t j = 0;
    while (i < dict.saD.size())
    {
      size_t left = i;

      auto sn = dict.saD[i];
      // Check if the suffix has length at least w and is not the complete phrase.
      auto phrase = dict.daD[i] + 1;
      assert(phrase > 0 && phrase < freq.size()); // + 1 because daD is 0-based
      size_t suffix_length = dict.select_b_d(dict.rank_b_d(sn + 1) + 1) - sn - 1;
      if (dict.b_d[sn] || suffix_length < w)
      {
        ++i; // Skip
      }
      else
      {
        // use the RMQ data structure to find how many of the following suffixes are the same except for the terminator (so they're the same suffix but in different phrases)
        // use the document array and the table of phrase frequencies to find the phrases frequencies and sum them up
        b_bwt[j++] = true;
        j += freq[phrase] - 1; // the next bits are 0s
        i++;
        if (i < dict.saD.size())
        {
          auto new_sn = dict.saD[i];
          auto new_phrase = dict.daD[i] + 1;
          assert(new_phrase > 0 && new_phrase < freq.size()); // + 1 because daD is 0-based
          size_t new_suffix_length = dict.select_b_d(dict.rank_b_d(new_sn + 1) + 1) - new_sn - 1;

          while (i < dict.saD.size() && (dict.lcpD[i] >= suffix_length) && (suffix_length == new_suffix_length))
          {
            j += freq[new_phrase];
            ++i;

            if (i < dict.saD.size())
            {
              new_sn = dict.saD[i];
              new_phrase = dict.daD[i] + 1;
              assert(new_phrase > 0 && new_phrase < freq.size()); // + 1 because daD is 0-based
              new_suffix_length = dict.select_b_d(dict.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
            }
          }
        }

        // Computing M
        size_t right = i - 1;
        M_entry_t m;
        m.len = suffix_length;
        m.left = dict.colex_daD[dict.rmq_colex_daD(left, right)];
        m.right = dict.colex_daD[dict.rMq_colex_daD(left, right)];

        M.push_back(m);
      }
    }
  }

};
#endif /* end of include guard: _PFP_HH */
