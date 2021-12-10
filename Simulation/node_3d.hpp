#ifndef NODE_3D_HPP
#define NODE_3D_HPP

#include "node_base.hpp"

struct node_3d : public node_base {

  bool within_limits(int x, int y, int z, int limit) {
    if(x < 0) return false;
    if(y < 0) return false;
    if(z < 0) return false;
    if(x >= limit) return false;
    if(y >= limit) return false;
    if(z >= limit) return false;

    return true;
  }


  void update_neighbors(std::vector< node_3d >& world,
                        size_t world_size) {

    static int relative_points[6][3] = { {-1,  0,  0},
                                         { 1,  0,  0},
                                         { 0,  1,  0},
                                         { 0, -1,  0},
                                         { 0,  0,  1},
                                         { 0,  0, -1} };

    for(int i = 0; i < 6; ++i) {
        int other_x = static_cast<int>(x_) + relative_points[i][0];
        int other_y = static_cast<int>(y_) + relative_points[i][1];
        int other_z = static_cast<int>(z_) + relative_points[i][2];
        if(within_limits(other_x, other_y, other_z, static_cast<int>(world_size))) {
            int other_pos = other_y +
                            (other_x + other_z * static_cast<int>(world_size))
                            * static_cast<int>(world_size);

            if(other_pos >= 0 && other_pos < static_cast<int>(world.size())) {
                node_3d* neighbor = &world[static_cast<size_t>(other_pos)];

                neighbors.push_back(neighbor);
              }
          }
      }
    inv_num_neighbors = 1.f / neighbors.size();
  }

  void set_coordinates(size_t row_size) {
    z_ = static_cast<float>(pos) / (row_size * row_size);
    size_t pos2 = pos - static_cast<size_t>(z_) * row_size * row_size;
    x_ = pos2 / row_size;
    y_ = pos2 % row_size;
  }
};



#endif // NODE_3D_HPP
