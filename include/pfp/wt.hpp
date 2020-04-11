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

class pfp_wt {
public:
  using wt_bv = sdsl::bit_vector;
  struct wt_node {
    wt_node (wt_node * parent = nullptr)
      : parent(parent)
    { }

    wt_bv bit_vector;
    wt_bv::rank_1_type bv_rank_1;
    wt_bv::select_1_type bv_select_1;
    wt_bv::select_0_type bv_select_0;

    uint32_t phrase_id = 0;

    wt_node * parent = nullptr;
    std::unique_ptr<wt_node> left;
    std::unique_ptr<wt_node> right;

    bool is_leaf () const {
      return !left.get();
    }

    bool is_root () const {
      return parent == nullptr;
    }
  };

  struct leaf_info {
    size_t alphabet_index;
    wt_node * leaf_link = nullptr;
  };

  pfp_wt()
    : root(new wt_node())
  { }

  pfp_wt(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse)
    : root(new wt_node()), alphabet(sorted_alphabet) {
    for (size_t i = 0; i < alphabet.size(); ++i) {
      leafs[alphabet[i]].alphabet_index = i;
    }

    create_bwt_rec(*root, alphabet.size(), 0, parse);
  }

  void construct(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse) {
    alphabet = sorted_alphabet;
    for (size_t i = 0; i < alphabet.size(); ++i) {
      leafs[alphabet[i]].alphabet_index = i;
    }

    create_bwt_rec(*root, alphabet.size(), 0, parse);
  }

  // public interface
  void print_leafs() {
    print_leafs_rec(root);
    std::cout << std::endl;
  }

  // operator[] - get phrase ID on i-th position
  uint32_t operator[] (const size_t i) {
    assert(i < root->bit_vector.size());

    return get_phrase_id(*root, i); // WT[] is 0-based
  }

  size_t rank(const size_t i, const uint32_t c) const {
    assert(i > 0 && i <= root->bit_vector.size());

    if (leafs.count(c) <= 0) {
      return 0;
    }

    wt_node * t_node = std::addressof(*root);
    size_t j = i;
    uint32_t alphabet_size = alphabet.size();
    uint32_t alphabet_start = 0;
    while (!t_node->is_leaf()) {
      const auto idx = leafs.at(c).alphabet_index - alphabet_start;
      const auto tres = alphabet_size / 2;
      const auto rank_bit_1 = t_node->bv_rank_1(j);

      if (idx < tres) {
        // left child
        j -= rank_bit_1;
        alphabet_size = tres;
        t_node = std::addressof(*t_node->left);
      }
      else {
        // right child
        j = rank_bit_1;
        alphabet_size -= tres;
        alphabet_start += tres;
        t_node = std::addressof(*t_node->right);
      }
    }

    return j;
  }

  size_t select(const size_t i, const uint32_t c) const {
    assert (i > 0 && i <= rank(size(), c));

    wt_node * t_node = leafs.at(c).leaf_link;
    size_t j = i;
    while (!t_node->is_root()) {
      if (std::addressof(*(t_node->parent->left)) == t_node) {
        j = t_node->parent->bv_select_0(j) + 1;
      }
      else {
        j = t_node->parent->bv_select_1(j) + 1;
      }
      t_node = t_node->parent;
    }

    return j - 1; // return index 0-based
  }

  size_t size() const {
    return root->bit_vector.size();
  }

  size_t range_count(const size_t t, const size_t b, const size_t i) {
    assert(i > 0 && i <= root->bit_vector.size());

    // t and b intervals of leafs
    // for each leaf from (t, b) - rank(i, c)
    size_t count = 0;
    for (size_t j = t; j <= b; ++j) {
      count += rank(i, alphabet[j]);
    }

    return count;
  }

  size_t range_select(const size_t t, const size_t b, const size_t r) {
    size_t lo = 1;
    size_t hi = size();
    assert (r > 0 && r <= range_count(t, b, hi));

    size_t ans = 0;

    // find r-th "filled column" in the interval (t, b)
    // binary search using rank_count
    while (lo <= hi) {
      const size_t i = (lo + hi) >> 1;
      const size_t m = range_count(t, b, i);

      if (m == r) {
        ans = i - 1;
        hi = i - 1;
      }

      if (m < r)
        lo = i + 1;
      else
        hi = i - 1;
    }

    return ans;
  }

private:
  // root node of wavelet tree
  std::unique_ptr<wt_node> root;
  std::vector<uint32_t> alphabet;
  std::map<uint32_t, leaf_info> leafs;

  uint32_t get_phrase_id(const wt_node & node, const size_t i) {
    if (node.is_leaf()) {
      return node.phrase_id;
    }

    const auto rank_bit_1 = node.bv_rank_1(i + 1);
    if (node.bit_vector[i]) {
      return get_phrase_id(*node.right, rank_bit_1 - 1);
    }
    else {
      return get_phrase_id(*node.left, i - rank_bit_1);
    }
  }

  void create_bwt_rec(wt_node & w_structure, const uint32_t alpha_size, const uint32_t alpha_start, const std::vector<uint32_t> & parse) {
    // divide alphabet to two sets
    const auto tres = alpha_size / 2;
    if (tres == 0) {
      // leaf
      w_structure.phrase_id = alphabet[alpha_start];
      leafs[w_structure.phrase_id].leaf_link = std::addressof(w_structure);

      return;
    }

    // create bit_vector
    w_structure.bit_vector.resize(parse.size());
    std::vector<uint32_t> parse_left;
    std::vector<uint32_t> parse_right;

    for (size_t i = 0; i < parse.size(); i++) {
      // TODO: direct access data structure
      const auto idx = leafs.at(parse[i]).alphabet_index - alpha_start;

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
    w_structure.bv_rank_1 = wt_bv::rank_1_type(&w_structure.bit_vector);
    w_structure.bv_select_1 = wt_bv::select_1_type(&w_structure.bit_vector);
    w_structure.bv_select_0 = wt_bv::select_0_type(&w_structure.bit_vector);

    w_structure.left = std::unique_ptr<wt_node>(new wt_node(std::addressof(w_structure)));
    w_structure.right = std::unique_ptr<wt_node>(new wt_node(std::addressof(w_structure)));

    create_bwt_rec(*(w_structure.left), tres, alpha_start, parse_left);
    create_bwt_rec(*(w_structure.right), alpha_size - tres, alpha_start + tres, parse_right);
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
