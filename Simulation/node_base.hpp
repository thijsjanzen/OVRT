#ifndef NODE_BASE_HPP
#define NODE_BASE_HPP

#include <vector>
#include <array>
#include <cmath>

enum cell_type {normal, cancer, infected, resistant, empty, max_num};

struct node_base {
private:
    cell_type node_type;

public:

  size_t pos;
  float x_;
  float y_;
  float z_; // not used in 2d
  size_t check_identifier;
  float inv_num_neighbors;
  float prob_normal_infected;

  float t_cell_concentration;
  float t_cell_death_rate;

  std::vector< node_base* > neighbors;

  node_base(node_base&&) = delete;
  const node_base& operator=(node_base&&) = delete;
  node_base(const node_base&) = delete;
  const node_base& operator=(const node_base&) = delete;

  node_base() {
    node_type = empty;
    pos = 0;
    x_ = 0;
    y_ = 0;
    inv_num_neighbors = 0.f;
    prob_normal_infected = 0.f;
    check_identifier = 0;
    t_cell_concentration = 0.f;
    t_cell_death_rate = 0.f;
  }

  void set_coordinates(size_t row_size);

  void update_neighbors(std::vector< node_base >& world,
                        size_t world_size);

  template<typename NODE>
  void add_neighbor(std::vector< NODE >& world,
                          size_t other_pos) {

      if(!neighbors.empty()) {
          for(auto i : neighbors) {
              if(i->pos == other_pos) {
                  return; // if the neighbor is already added, don't add again.
              }
          }
      }

      NODE* neighbor = &world[static_cast<size_t>(other_pos)];
      neighbors.push_back(neighbor);
  }

  std::array<float, 4> calc_prob_of_growth() const {
    std::array<float, 4> prob_of_growth = {0.f, 0.f, 0.f, 0.f};

    if(node_type == normal) {
      // infected can grow into this, but at lower frequency:
      prob_of_growth[infected]  =  prob_normal_infected * freq_type_neighbours(infected);
    }
    if(node_type == cancer) {
      // infected nodes can grow into this
      prob_of_growth[infected]  = freq_type_neighbours(infected);
    }
    if(node_type == empty) {
      prob_of_growth[normal]    = freq_type_neighbours(normal);
      prob_of_growth[cancer]    = freq_type_neighbours(cancer);
      prob_of_growth[resistant] = freq_type_neighbours(resistant);
    }

    return prob_of_growth;
  }

  float freq_type_neighbours(const cell_type& ref_type) const {

    auto count = std::count_if(neighbors.begin(), neighbors.end(),
                  [&ref_type](auto i){return ref_type == i->get_cell_type(); });

    return static_cast<float>(count * inv_num_neighbors);
  }

  std::vector< size_t > get_cancer_neighbours() const {
    std::vector< size_t > output;
    for(auto i : neighbors) {
      if(i->get_cell_type() == cancer) {
        output.push_back(i->pos);
      }
    }
    return output;
  }

  cell_type get_cell_type() const {
    return node_type;
  }

  void set_cell_type(cell_type new_type) {
    node_type = new_type;
  }

  size_t get_identifier() {
    return check_identifier;
  }

  void set_identifier(size_t id) {
    check_identifier = id;
  }

  void add_t_cell(float amount) {
    t_cell_concentration += amount;
  }

  void update_t_cell_death_rate(float t_cell_rate,
                                float t_cell_density_scaler,
                                float t_cell_inflection_point) {
    if (t_cell_concentration < 1e-5f) {
        t_cell_death_rate = 0.f;
        return;
     }

    float denominator = 1.f + expf(5 * (t_cell_inflection_point - t_cell_concentration)); // same as -1 * (b - mu)
    float added_t_cell_death_rate = t_cell_rate * 1.f / denominator;

    if (added_t_cell_death_rate >= std:: numeric_limits<float>::max()) {
        added_t_cell_death_rate = 1e9;
     }

    float mult = 1.0f - t_cell_density_scaler *
        freq_type_neighbours(cancer);

    if(mult < 0.f) mult = 0.f;
    t_cell_death_rate = mult * added_t_cell_death_rate;
    if (std::isinf(t_cell_death_rate)) {
        t_cell_death_rate = 1e9;
      }
    if (t_cell_death_rate >= std:: numeric_limits<float>::max()) {
        t_cell_death_rate = 1e9;
      }
    return;
  }

};

#endif // NODE_BASE_HPP
