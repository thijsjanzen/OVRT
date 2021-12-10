#ifndef VORONOI_TOOLS_HPP
#define VORONOI_TOOLS_HPP
//
//  voronoi_tools.hpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#include <cstdio>
#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include "voronoi.hpp"

struct voronoi_point {
    double x_, y_;

    voronoi_point() {
        x_ = 0.0; y_ = 0.0;
    }

    voronoi_point& operator=(const voronoi_point& other) {
      x_ = other.x_;
      y_ = other.y_;
      return *this;
    }

    voronoi_point(const voronoi_point& other) {
        x_ = other.x_;
        y_ = other.y_;
    }

    voronoi_point(double x, double y) : x_(x), y_(y) {}

    bool operator==(const voronoi_point& other) const {
        // explicitly, this doesn't keep track of left and right!
      //  if(x_ != other.x_) return false;
      //  if(y_ != other.y_) return false;
      if(abs(x_ - other.x_) > 1e-3) return false;
      if(abs(y_ - other.y_) > 1e-3) return false;

      return true;
    }

    bool operator<(const voronoi_point& other) const {
       // if(fabs(x_ - other.x_) < 1e-4f) return y_ < other.y_;
       // return x_ < other.x_;
       if(x_ == other.x_) return y_ < other.y_;
       return x_ < other.x_;
    }

    bool operator!=(const voronoi_point& other) const {
        return !(*this == other);
    }
};

struct voronoi_edge {
  voronoi_edge() {
  }

    voronoi_edge(voronoi_point s, voronoi_point e,
                 size_t l, size_t r) : start(s), end(e), left(l), right(r) {}

    voronoi_point start;
    voronoi_point end;
    size_t left;
    size_t right;

    voronoi_edge& operator=(const voronoi_edge& other) {
      left  = other.left;
      right = other.right;
      start = other.start;
      end   = other.end;
      return *this;
    }

    voronoi_edge(const voronoi_edge& other) {
      left  = other.left;
      right = other.right;
      start = other.start;
      end   = other.end;
    }

    bool operator<(const voronoi_edge& other) const {
      return start < other.start;
    }

    bool operator==(const voronoi_edge& other) const {
     if(start != other.start) return false;
     if(end   != other.end)   return false;
     return true;
    }

    bool operator!=(const voronoi_edge& other) const {
        return !(*this == other);
    }

    bool check() const {
      if(std::isnan(start.x_)) {
          std::cout << "start.x_ isnan\n";
          return false;
      }
      if(std::isnan(start.y_)) {
          std::cout << "start.y_ isnan\n";
          return false;
      }
      if(std::isnan(end.x_)) {
          std::cout << "end.x_ isnan\n";
          return false;
      }
      if(std::isnan(end.y_)) {
          std::cout << "end.y_ isnan\n";
          return false;
      }
      return true;
    }

    double calc_dist() {
      double x_x = (start.x_ - end.x_) * (start.x_ - end.x_);
      double y_y = (start.y_ - end.y_) * (start.y_ - end.y_);
      return(std::sqrt(x_x + y_y));
    }
};

namespace voronoi_tools {

  inline std::vector< double > calc_dist(const std::vector< voronoi_edge >& v,
                                  const voronoi_edge& point) {
    std::vector< double > output(v.size());
    for(size_t i = 0; i < v.size(); ++i) {
      double dx = point.end.x_ - v[i].start.x_;
      double dy = point.end.y_ - v[i].start.y_;
      output[i] = dx * dx + dy * dy;
    }
    return output;
  }

  inline void invert_edges(std::vector< voronoi_edge>& edges, size_t pos) {
      // check for inverted edges.
      for(auto& i : edges) {
          if(i.right != pos) {
              size_t tmp = i.right;
              i.right = i.left;
              i.left = tmp;

              auto tmp_2 = i.start;
              i.start = i.end;
              i.end = tmp_2;
          }
      }
  }

  inline std::vector< voronoi_point> clean_edges(const std::vector< voronoi_edge >& input_edges,
                                          size_t pos) {
      // the goal is to connect all edges, and then provide
      // all the outer points
      std::vector< voronoi_edge > edges = input_edges;
      invert_edges(edges, pos);

      std::vector< voronoi_edge > new_edges;

      voronoi_edge focal_edge = edges.back();
      new_edges.push_back(focal_edge);

      while(!edges.empty()) {
         std::vector< double > distances = calc_dist(edges, focal_edge);
         auto min_val = std::min_element(distances.begin(), distances.end());
         auto match = std::distance(distances.begin(), min_val);

         focal_edge = edges[match];
         edges[match] = edges.back();
         edges.pop_back();
         new_edges.push_back(focal_edge);
      }
      // allright, we have connected them now
      // now we collect the starting points
      std::vector< voronoi_point > outer_points;
      for(auto i : new_edges) {
          voronoi_point temp_point = i.start;
          outer_points.push_back(temp_point);
      }
      return outer_points;
  }
}


#endif // VORONOI_TOOLS_HPP
