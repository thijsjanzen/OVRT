#ifndef RANDOM_THIJS_H
#define RANDOM_THIJS_H

// Copyright 2019 Thijs Janzen

#include <random>
#include <map>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include "rndutils.hpp"

struct rnd_t {
  rndutils::xorshift128 rndgen_;

  rnd_t() {
    std::random_device rd;
    rndgen_ = rndutils::make_random_engine(rd());
  }

  rndutils::uniform01_distribution<float> rndutil_norm;

  size_t random_number(size_t n)    {
    if(n <= 1) return 0;
    return std::uniform_int_distribution<size_t> (0, n - 1)(rndgen_);
  }

  void set_seed(size_t s) {
    rndgen_ = rndutils::make_random_engine(s);
  }

  float uniform()    {
    return rndutil_norm(rndgen_);
  }

  float Expon(float lambda) {
    if (lambda == 0.f) return 1e20f;
    return std::exponential_distribution<float>(lambda)(rndgen_);
  }
};

struct binned_distribution {
public:
  binned_distribution() {
  }

  binned_distribution(size_t num_of_bins) : num_bins(num_of_bins) {
    bin_size = num_of_bins;
    row_sum.resize(num_bins);
  }

  binned_distribution(size_t num_of_bins, size_t num_values) {
    bin_size = num_values / num_of_bins;
    num_bins = num_of_bins;

    row_sum.resize(num_bins, 0.f);
    values.resize(num_values, 0.f);
  }

  template< typename It>
  binned_distribution(It first, It last, size_t num_of_bins) :
  num_bins(num_of_bins) {
    size_t N = std::distance(first, last);
    bin_size = N / num_bins;
    row_sum.resize(num_bins);
    values = std::vector<float>(first, last);

    for(size_t i = 0; i < num_bins; ++i) {
      auto start_it = values.begin() + i * bin_size;
      auto end_it = start_it + bin_size;
      row_sum[i] = std::accumulate(start_it, end_it, 0.f);
    }
    total_sum = std::accumulate(row_sum.begin(), row_sum.end(), 0.f);
  }

  template <typename It>
  size_t draw_from_dist(It first, It last, float max_val, rnd_t& r) const {
   size_t max_index = static_cast<size_t>(std::distance(first, last));
   size_t cnt = 0;

   if (max_val == 0.f) {
      return r.random_number(static_cast<size_t>(max_index));
   }

   while (true) {
     size_t index = r.random_number(static_cast<size_t>(max_index));
     float val = *(first + index);
     if (r.uniform() < (1.f * val / max_val)) {
       return index;
     }
     cnt++;
     if (cnt > max_index * 10) {
       return draw_cdf(first, last, r);
     }
   }
  }

  size_t draw_explicit(rnd_t& r) const {
    size_t row = draw_cdf(row_sum.begin(), row_sum.end(), r);
    size_t col;
    float frac = row_sum[row] * 1.f / bin_size;
    auto start = values.begin() + row * bin_size;
    auto end = start + bin_size;

    if(frac < 0.1f) {
      col = draw_cdf(start, end, r);
    } else {
      auto max_val = 1.0; //std::max_element(start, end);
      col = draw_from_dist(start, end, max_val, r);
    }
    size_t result = row * bin_size + col;

    return result;
  }

  template< typename It>
  size_t draw_cdf(It first, It last, rnd_t& r) const {
    auto draw_dist = rndutils::mutable_discrete_distribution< int, rndutils::all_zero_policy_nothing >{};
    draw_dist.mutate(first, last);
    return(static_cast<size_t>(draw_dist(r.rndgen_)));
  }

  void update_entry(size_t pos, float new_val) {
    float old_val = values[pos];

    if(old_val == new_val) return;
    values[pos] = new_val;
    size_t row = pos / bin_size;

    if (new_val - old_val > 1e6f) {
      double diff = static_cast<double>(new_val) - static_cast<double>(old_val);
      row_sum[row] += static_cast<float>(diff);
      total_sum += static_cast<float>(diff);
    } else {
        float diff = new_val - old_val;
        row_sum[row] += diff;
        total_sum += diff;
     }


    if (total_sum < 0.f) total_sum = 0.f;
    if (total_sum > 1e9f) total_sum = 1e9f;

    if(row_sum[row] < 0.f) row_sum[row] = 0.f;
    if (row_sum[row] > 1e9f) row_sum[row] = 1e9f;

  }

  float get_total_sum() const {
    return total_sum;
  }

  void update_all() {
    for(size_t row = 0; row < num_bins; ++row) {
      auto start = values.begin() + row * bin_size;
      auto end = start + bin_size;
      row_sum[row] = std::accumulate(start, end, 0.f);
      if(row_sum[row] < 0.f) row_sum[row] = 0.f;
    }
    total_sum = std::accumulate(row_sum.begin(), row_sum.end(), 0.f);
  }

  float get_value(size_t pos) const {
    return values[pos];
  }

private:
  float total_sum;
  size_t bin_size;
  size_t num_bins;
  std::vector<float> row_sum;
  std::vector<float> values;
};

#endif  // RANDOM_THIJS_H
