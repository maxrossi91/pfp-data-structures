/* pfp-parse - prefix free parsing parse
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
   \file parse.hpp
   \brief parse.hpp define and build the prefix-free parse data structure.
   \author Massimiliano Rossi
   \date 03/04/2020
*/


#ifndef _PFP_PARSE_HH
#define _PFP_PARSE_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C" {
    #include<gsacak.h>
}

// TODO: Extend it to non-integer alphabets
class parse{
public:
  std::vector<uint32_t> p;
  std::vector<uint_t> saP;
  std::vector<uint_t> isaP;
  std::vector<int_t> lcpP;
  sdsl::rmq_succinct_sct<> rmq_lcp_P;
  // sdsl::bit_vector b_p; // Starting position of each phrase in D
  // sdsl::bit_vector::rank_1_type rank_b_p;
  // sdsl::bit_vector::select_1_type select_b_p;
  bool saP_flag = false;
  bool isaP_flag = false;
  bool lcpP_flag = false;
  bool rmq_lcp_P_flag = false;

  size_t alphabet_size;

  typedef size_t size_type;

  parse(  std::vector<uint32_t>& p_,
          size_t alphabet_size_,
          bool saP_flag_ = true,
          bool isaP_flag_ = true,
          bool lcpP_flag_ = true,
          bool rmq_lcp_P_flag_ = true ):
          p(p_),
          alphabet_size(alphabet_size_)
  {
    assert(p.back() == 0);
    build(saP_flag_, isaP_flag_, lcpP_flag_, rmq_lcp_P_flag_);


  }

  parse(  std::string filename,
          size_t alphabet_size_,
          bool saP_flag_ = true,
          bool isaP_flag_ = true,
          bool lcpP_flag_ = true,
          bool rmq_lcp_P_flag_ = true ):
          alphabet_size(alphabet_size_)
  {
    // Building dictionary from file
    std::string tmp_filename = filename + std::string(".parse");
    read_file(tmp_filename.c_str(), p);
    p.push_back(0); // this is the terminator for the sacak algorithm

    build(saP_flag_, isaP_flag_, lcpP_flag_, rmq_lcp_P_flag_);

  }

  void build(bool saP_flag_, bool isaP_flag_, bool lcpP_flag_, bool rmq_lcp_P_flag_){

    // TODO: check if it has been already computed
    if(saP_flag_){
      saP.resize(p.size());
      // suffix array of the parsing.
      verbose("Computing SA of the parsing");
      _elapsed_time(
        sacak_int(&p[0],&saP[0],p.size(),alphabet_size);
      );
    }

    assert(!isaP_flag_ || (saP_flag || saP_flag_) );
    if(isaP_flag_ && !isaP_flag){
      // inverse suffix array of the parsing.
      verbose("Computing ISA of the parsing");
      _elapsed_time(
        {
          isaP.resize(p.size());
          for(int i = 0; i < saP.size(); ++i){
            isaP[saP[i]] = i;
          }
        }
      );

    }

    if(lcpP_flag_){
      lcpP.resize(p.size());
      // LCP array of the parsing.
      verbose("Computing LCP of the parsing");
      _elapsed_time(
        LCP_array(&p[0], isaP, saP, p.size(), lcpP);
      );
    }


    assert(!rmq_lcp_P_flag_ || (lcpP_flag || lcpP_flag_));
    if(rmq_lcp_P_flag_ && ! rmq_lcp_P_flag){
      rmq_lcp_P_flag = true;
      verbose("Computing RMQ over LCP of the parsing");
      // Compute the LCP rank of P
      _elapsed_time(
        rmq_lcp_P = sdsl::rmq_succinct_sct<>(&lcpP);
      );
    }

  }

  // Serialize to a stream.
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += sdsl::serialize(p, out, child, "parse");
    written_bytes += sdsl::serialize(saP, out, child, "saP");
    written_bytes += sdsl::serialize(isaP, out, child, "isaP");
    written_bytes += sdsl::serialize(lcpP, out, child, "lcpP");
    written_bytes += rmq_lcp_P.serialize(out, child, "rmq_lcp_P");
    // written_bytes += b_p.serialize(out, child, "b_p");
    // written_bytes += rank_b_p.serialize(out, child, "rank_b_p");
    // written_bytes += select_b_p.serialize(out, child, "select_b_p");
    written_bytes += sdsl::write_member(alphabet_size, out, child, "alphabet_size");

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  //! Load from a stream.
  void load(std::istream &in)
  {
    sdsl::load(p, in);
    sdsl::load(saP, in);
    sdsl::load(isaP, in);
    sdsl::load(lcpP, in);
    rmq_lcp_P.load(in);
    // b_p.load(in);
    // rank_b_p.load(in);
    // select_b_p.load(in);
    sdsl::read_member(alphabet_size, in);
  }

};

#endif /* end of include guard: _PFP_PARSE_HH */
