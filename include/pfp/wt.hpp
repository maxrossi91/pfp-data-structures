/* pfp - prefix free parsing wt W support
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
   \file wt.hpp
   \brief wt.hpp define and build the prefix-free parsing wavelet tree support data structures.
   \author Ondřej Cvacho
   \date 03/04/2020
*/


#ifndef _PFP_WT_HH
#define _PFP_WT_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include<pfp.hpp>

class pfp_wt {
public:
  using wt_bv = sdsl::bit_vector;
  struct wt_node {
    wt_bv bit_vector;
    wt_bv::rank_1_type bv_rank;
    wt_bv::select_1_type bv_select;

    uint32_t phrase_id;

    std::unique_ptr<wt_node> left;
    std::unique_ptr<wt_node> right;

    bool is_leaf () {
      return !left.get();
    }
  };

  pfp_wt()
    : root(new wt_node())
  { }

  pfp_wt(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse)
    : root(new wt_node()) {
    verbose("Construction of WT of P");
    create_bwt_rec(*root, sorted_alphabet, parse);
  }

  void construct(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse) {
    verbose("Construction of WT of P");
    create_bwt_rec(*root, sorted_alphabet, parse);
  }

  // public interface
  void print_leafs() {
    print_leafs_rec(root);
    std::cout << std::endl;
  }

  // rank, select

private:
  // root node of wavelet tree
  std::unique_ptr<wt_node> root;

  void create_bwt_rec(wt_node & w_structure, const std::vector<uint32_t> & alphabet, const std::vector<uint32_t> & parse) {
    // divide alphabet to two sets
    const uint32_t tres = alphabet.size() / 2;
    if (tres == 0) {
      // leaf
      w_structure.phrase_id = alphabet[0];

      return;
    }

    // create bit_vector
    w_structure.bit_vector.resize(parse.size());
    std::vector<uint32_t> parse_left;
    std::vector<uint32_t> parse_right;

    // std::cout << "> Creating bit vector" << std::endl;
    for (size_t i = 0; i < parse.size(); i++) {
      // TODO: direct access data structure
      const auto iter = std::find(alphabet.begin(), alphabet.end(), parse[i]);
      const auto idx = std::distance(alphabet.begin(), iter);

      if (idx < tres) {
        w_structure.bit_vector[i] = 0;
        parse_left.push_back(parse[i]);
      }
      else {
        w_structure.bit_vector[i] = 1;
        parse_right.push_back(parse[i]);
      }
    }

    // rank & select support
    w_structure.bv_rank = wt_bv::rank_1_type(&w_structure.bit_vector);
    w_structure.bv_select = wt_bv::select_1_type(&w_structure.bit_vector);

    w_structure.left = std::unique_ptr<wt_node>(new wt_node());
    w_structure.right = std::unique_ptr<wt_node>(new wt_node());

    create_bwt_rec(*(w_structure.left), std::vector<uint32_t>(alphabet.begin(), alphabet.begin() + tres), parse_left);
    create_bwt_rec(*(w_structure.right), std::vector<uint32_t>(alphabet.begin() + tres, alphabet.end()), parse_right);
  }

  void print_leafs_rec(const std::unique_ptr<wt_node> & node) {
    if (node->is_leaf()) {
      std::cout << node->phrase_id << " | ";
    }
    else {
      print_leafs_rec(node->left);
      print_leafs_rec(node->right);
    }
  }
};

#endif /* end of include guard: _PFP_WT_HH */
