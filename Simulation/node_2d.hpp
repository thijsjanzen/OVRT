#ifndef NODE_2D_HPP
#define NODE_2D_HPP

#include "node_base.hpp"

struct node_2d : public node_base {

  void update_neighbors(std::vector< node_2d >& world,
                              size_t world_size) {

      static int relative_points[4][2] = { {-1, 0},
                                          {1, 0},
                                          {0, 1},
                                          {0, -1} };

      for(int i = 0; i < 4; ++i) {
          int other_x = static_cast<int>(x_) + relative_points[i][0];
          int other_y = static_cast<int>(y_) + relative_points[i][1];
          if(other_x >= 0 &&
             other_y >= 0 &&
             other_x < static_cast<int>(world_size) &&
             other_y < static_cast<int>(world_size)) {
              int other_pos = other_y + other_x * static_cast<int>(world_size);
              if(other_pos >= 0 && other_pos < static_cast<int>(world.size())) {
                  node_2d* neighbor = &world[static_cast<size_t>(other_pos)];

                  neighbors.push_back(neighbor);
              }
          }
      }
      inv_num_neighbors = 1.f / neighbors.size();
  }

  void set_coordinates(size_t row_size) {
      x_ = pos / row_size;
      y_ = pos % row_size;
  }
};



#endif // NODE_2D_HPP
