//
//  simulation.hpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#ifndef simulation_hpp
#define simulation_hpp

#include <cstdio>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <memory>
#include <chrono>
#include <fstream>
#include <iostream>

#include "parameters.hpp"
#include "random_thijs.hpp"
#include "node_3d.hpp"
#include "node_2d.hpp"
#include "voronoi_tools.hpp"

class simulation {
public:
  float t;
  size_t num_cells;
  size_t sq_size;
  float total_t_cell_concentration;
  Param parameters;
  virtual ~simulation() {}
  virtual void update_one_step() = 0;

  virtual cell_type get_cell_type(size_t) const = 0;
  virtual float get_t_cell_concentration(size_t) const = 0;
  virtual int get_nni(size_t pos) const = 0;


  virtual std::array<size_t, 5> get_count_cell_types() const = 0;
  virtual void initialize_network(std::vector< std::vector< voronoi_point > >& all_polys,
                                  grid_type used_grid_type) = 0;

  virtual const std::array< binned_distribution, 4 >& get_growth_prob()  = 0;
  virtual const std::array< binned_distribution, 4 >& get_death_prob()  = 0;
  virtual const std::array< binned_distribution, 4 >& get_T_cell_death_prob() = 0;

  virtual void set_infection_type(infection_routine) = 0;
  virtual void set_percent_infected(float percent_infected) = 0;
  virtual void set_start_setup(start_type new_type) = 0;
  virtual void add_infected(infection_routine infect_type,
                                 float fraction) = 0;
  virtual void update_grow_params() = 0;
  virtual infection_routine get_infection_type() const = 0;
  virtual float get_percent_infected() const = 0;

  virtual void add_cells(const cell_type& focal_cell_type) = 0;

  virtual void obtain_equilibrium(bool verbose) = 0;

  virtual float get_t_cell_rate(size_t) = 0;
  virtual void reset_t_cell_death_rate() = 0;
};


template <typename NODE>
class simulation_impl : public simulation {
public:
  bool using_3d;

  std::vector< NODE > world;
  rnd_t rndgen;

  std::array< binned_distribution, 4 > growth_prob;
  std::array< binned_distribution, 4 > death_prob;
  std::array< binned_distribution, 4 > t_cell_death_prob;

  const std::array< binned_distribution, 4 >& get_growth_prob() override {
    return growth_prob;
  }

  const std::array< binned_distribution, 4 >& get_death_prob()  override {
      return death_prob;
  }

  const std::array< binned_distribution, 4 >& get_T_cell_death_prob()  override {
      return t_cell_death_prob;
  }

  std::array<size_t, 5> num_cell_types;

  simulation_impl(const Param& param) :
    world(param.sq_num_cells * param.sq_num_cells)
  {
    parameters = param;
    using_3d = false;
    rndgen.set_seed(parameters.seed);
    total_t_cell_concentration = 0.f;
    sq_size = parameters.sq_num_cells;
    num_cells = sq_size * sq_size;

    num_cell_types = {0, 0, 0, 0, num_cells}; // all cells are empty

    for (size_t i = 0; i < world.size(); ++i) {
        world[i].pos = i;
        world[i].set_coordinates(sq_size);
        world[i].prob_normal_infected = parameters.prob_normal_infection;
      }

    binned_distribution temp(sq_size, num_cells);
    for(size_t i = 0; i < 4; ++i) {
        growth_prob[i] = temp;
        death_prob[i] = temp;
        t_cell_death_prob[i] = temp;
      }

    long_distance_infection_probability = std::vector<double>(sq_size, parameters.prob_infection_upon_death);
  }

  simulation_impl(const Param& param,
             bool have_to_use_3d) :
    world(param.sq_num_cells * param.sq_num_cells * param.sq_num_cells)
  {
    parameters = param;
    using_3d = have_to_use_3d;
    rndgen.set_seed(parameters.seed);

    sq_size = parameters.sq_num_cells;
    num_cells = sq_size * sq_size * sq_size;\

    num_cell_types = {0, 0, 0, 0, num_cells}; // all cells are empty

    for (size_t i = 0; i < world.size(); ++i) {
        world[i].pos = i;
        world[i].set_coordinates(sq_size);
        world[i].prob_normal_infected = parameters.prob_normal_infection;
      }

    binned_distribution temp(sq_size, num_cells);
    for(size_t i = 0; i < 4; ++i) {
        growth_prob[i] = temp;
        death_prob[i] = temp;
        t_cell_death_prob[i] = temp;
      }

    long_distance_infection_probability = std::vector<double>(sq_size, parameters.prob_infection_upon_death);
  }

  void update_one_step() override {
    update_rates();
    float lambda = std::accumulate(rates.begin(), rates.end(), 0.0f);

    float dt = rndgen.Expon(lambda);
    while(std::isinf(dt)) {
        static int counter = 0;
        dt = rndgen.Expon(lambda);
        counter++;
        if(counter > 100) {
            exit(1);
          }
      }

    size_t event = pick_event(rates, lambda);
    do_event(event);

    if(parameters.start_setup == grow) {
        if(t < parameters.time_adding_cancer &&
           t+dt >= parameters.time_adding_cancer) {
            add_cells(cancer);
          }

        if(t < parameters.time_adding_virus &&
           t+dt >= parameters.time_adding_virus) {
            add_infected(parameters.infection_type,
                         parameters.percent_infected);
          }

        if(t < parameters.time_adding_virus_2 &&
           t+dt >= parameters.time_adding_virus_2) {
            add_infected(parameters.infection_type_2,
                         parameters.percent_infected_2);
          }
      }

    // check if the distributions are not getting close to zero,
    // in that case, some numerical irregularities might pop up
    // check only every hour:
    int delta_t = static_cast<int>(t + dt) - static_cast<int>(t);
    if (delta_t > 0) {
        for(size_t i = 0; i < 4; ++i) {
            if (num_cell_types[i] < 100 &&
                growth_prob[i].get_total_sum() != 0.f) {
                growth_prob[i].update_all();
              }
          }
      }

    if(total_t_cell_concentration > 0.f) {
        int delta_t_ms = static_cast<int>(10 * (t + dt)) - static_cast<int>(10 * t);
        if(delta_t_ms > 0) {
            diffuse();
        }

        /*if(delta_t_ms > 0) {
           reset_nni();
           update_nni();
        }*/
      }

    t += dt;
  }

  /*void reset_nni() {
      for (auto& i : world) {
          i.nearest_cancer_cell = -1;
      }
  }

  void update_nni() {
      for (auto& i : world) {
          if (i.get_cell_type() == cell_type::infected) {
              auto nni = i.find_nearest_cancer_cell(rndgen.random_number(1e9));
              i.nearest_cancer_cell = nni;
          }
      }
  }*/

  void diffuse() {
    std::vector<float> new_concentration(world.size(), 0.f);
    for(size_t i = 0; i < world.size(); ++i) {
      if(world[i].t_cell_concentration > 0.f) {

        float current_conc = world[i].t_cell_concentration;
        new_concentration[i] += current_conc;

        for(const auto& j : world[i].neighbors) {
           size_t other_pos = j->pos;
           float other_conc = world[other_pos].t_cell_concentration;
           float delta_conc = current_conc - other_conc;

           float diffusion_amount = delta_conc *
               parameters.diffusion * world[i].inv_num_neighbors;

           if (diffusion_amount > 0) { // otherwise we track the same flow twice.
             new_concentration[i] -= diffusion_amount;

             new_concentration[other_pos] += diffusion_amount;
           }
        }
      }
    }

    for(size_t i = 0; i < new_concentration.size(); ++i) {
        if (new_concentration[i] > 0.f) {
            new_concentration[i] *= (1.f - parameters.evaporation);
            if(new_concentration[i] < 1e-5f) new_concentration[i] = 0.f;
        }
        world[i].t_cell_concentration = new_concentration[i];  // swap of the vectors
        world[i].update_t_cell_death_rate(parameters.t_cell_rate,
                                          parameters.t_cell_density_scaler,
                                          parameters.t_cell_inflection_point);

        update_t_cell_death_prob(world[i].t_cell_death_rate, i);
    }
    total_t_cell_concentration = std::accumulate(new_concentration.begin(),
                                                 new_concentration.end(),
                                                 0.f);
  }

  void increase_t_cell_concentration(size_t pos) {
    world[pos].add_t_cell(parameters.t_cell_increase);
    total_t_cell_concentration += parameters.t_cell_increase;
  }

  cell_type get_cell_type(size_t pos) const override {
    return world[pos].get_cell_type();
  }

  float get_t_cell_concentration(size_t pos) const override {
    return world[pos].t_cell_concentration;
  }

  int get_nni(size_t pos) const override {
      return world[pos].nearest_cancer_cell;
  }

  void update_growth_prob(size_t pos) {
    std::array<float, 4> probs = world[pos].calc_prob_of_growth();
    for (size_t i = 0; i < 4; ++i) {
        growth_prob[i].update_entry(pos, probs[i]);
      }
  }

  void update_death_prob(size_t pos) {
    for (size_t i = 0; i < 4; ++i) {
        float new_val = 0.f;
        if (i == world[pos].get_cell_type()) new_val = 1.f;

        death_prob[i].update_entry(pos, new_val);
      }
  }

  void update_t_cell_death_prob(float new_prob, size_t pos) {
    for (size_t i = 0; i < 4; ++i) {
        float new_val = 0.f;
        if (i == world[pos].get_cell_type()) new_val = new_prob;

        t_cell_death_prob[i].update_entry(pos, new_val);
      }
  }


  void update_death_prob(size_t pos,
                         cell_type old_type,
                         cell_type new_type) {
    if(old_type < death_prob.size()) death_prob[old_type].update_entry(pos, 0.f);
    if(new_type < death_prob.size()) {
        float death_rate = 1.0f;

        death_prob[new_type].update_entry(pos, death_rate);
      }
  }


  void set_percent_infected(float percent_infected) override {
    parameters.percent_infected = percent_infected;
  }

  void set_infection_type(infection_routine infect_routine) override {
    parameters.infection_type = infect_routine;
  }

  std::array<size_t, 5> get_count_cell_types() const override {
    return num_cell_types;
  }

  void add_infected(infection_routine infect_type,
                    float fraction) override {
    if(fraction == 1.0f &&
       infect_type == random_infection) {
        infect_all_cancer();
        return;
      }

    switch(infect_type) {
      case random_infection:
        infect_random(fraction);
        break;
      case center_infection:
        infect_center_largest(fraction);
        break;
      case periphery_infection:
        infect_periphery(fraction);
        break;
      }
  }

  std::array<size_t, 5> count_cell_types() const {
    std::array<size_t, 5> total_num_cell_types = {0, 0, 0, 0, 0};
    for(const auto& i : world) {
        total_num_cell_types[ i.get_cell_type()]++;
      }
    return total_num_cell_types;
  }

  void add_cells(const cell_type& focal_cell_type) override {

    // pick center node:
    cell_type to_be_replaced = empty;
    if (focal_cell_type == cancer  ) to_be_replaced = normal;
    if (focal_cell_type == infected) to_be_replaced = cancer;

    if(num_cell_types[to_be_replaced] == 0) {
        return;
      }

    size_t focal_pos = find_central_cell(to_be_replaced);

    std::vector<size_t> cells_turned(1, focal_pos);

    change_cell_type(focal_pos, focal_cell_type);

    auto max_number_of_cells = static_cast<size_t>(parameters.initial_number_cancer_cells);

    if (focal_cell_type == normal) max_number_of_cells = parameters.initial_number_normal_cells;

    if (max_number_of_cells > world.size()) max_number_of_cells = world.size();

    if(parameters.start_setup == converge && focal_cell_type == normal) {
        max_number_of_cells = world.size() * 0.95f;
      }

    size_t counter = 0;
    while (cells_turned.size() < max_number_of_cells &&
           counter < cells_turned.size()) {
        focal_pos = cells_turned[counter];
        for (size_t i = 0; i < world[focal_pos].neighbors.size(); ++i) {
            size_t other_pos = world[focal_pos].neighbors[i]->pos;
            if (world[other_pos].get_cell_type() != focal_cell_type) {
                change_cell_type(other_pos, focal_cell_type);
                cells_turned.push_back(other_pos);
                if (cells_turned.size() >= max_number_of_cells) break;
              }
          }
        counter++;
      }

    for(size_t i = 0; i < num_cells; ++i) {
        update_growth_prob(i);
        update_death_prob(i);
      }
  }

  void obtain_equilibrium(bool verbose = false) override {

    std::vector< float > densities(10, 0);
    t = 0.f;
    float prev_t = t;
    std::array<size_t, 5> cell_counts = num_cell_types;
    int count = 0;
    while(t < parameters.maximum_time) {
        update_one_step();

        if(static_cast<int>(t) - static_cast<int>(prev_t) == 10) {
            cell_counts = num_cell_types;
            auto density_normal = 1.f * cell_counts[normal] / (sq_size * sq_size);
            densities[count % 10] = cell_counts[normal];
            count++;
            if(count / 10 > 1) {
                float sum_first_half = 0.f;
                float sum_second_half = 0.f;
                for(size_t i = 0; i < 5; ++i) {
                    sum_first_half  += densities[i ];
                    sum_second_half += densities[i + 5];
                  }

                if(sum_first_half >= sum_second_half && density_normal > 0.5) {
                    break;
                  }
              }
            prev_t = t;
            if (verbose) {
                std::cout << t << " " << cell_counts[normal] << "\n";
              }
          }
      }
    return;
  }

  Param get_parameters() {
    return parameters;
  }

  infection_routine get_infection_type() const override {
    return(parameters.infection_type);
  }
  float get_percent_infected() const override {
    return(parameters.percent_infected);
  }
  void set_start_setup(start_type new_type) override {
    parameters.start_setup = new_type;
  }

  float get_t_cell_rate(size_t pos) override {
      auto local_cell_type = world[pos].get_cell_type();
      if (local_cell_type == empty) return 0.f;

      auto local_death_rate = 0.0;

      if (parameters.t_cell_sensitivity[local_cell_type]) {
          local_death_rate = t_cell_death_prob[local_cell_type].get_value(pos);
      }

      return local_death_rate;
  }

  void reset_t_cell_death_rate() override {
    for (size_t pos = 0; pos < world.size(); ++pos) {
        world[pos].t_cell_concentration = 0.0;
        world[pos].t_cell_death_rate = 0.0;
        auto local_cell_type = world[pos].get_cell_type();
        if (local_cell_type != empty) {
            t_cell_death_prob[local_cell_type].update_entry(pos, 0.0);
        }
    }
  }

  // size_t find_center(const cell_type& focal_cell_type) const;
  size_t find_central_cell(const cell_type& focal_cell_type) const {
    // first calculate average x and y of cell type:

    float x = 0.f;
    float y = 0.f;
    float z = 0.f;
    int counter = 0;
    for(const auto& i : world) {
        if(i.get_cell_type() == focal_cell_type) {
            x += i.x_;
            y += i.y_;
            if (using_3d) z += i.z_;
            counter++;
          }
      }
    x *= 1.0f / counter;
    y *= 1.0f / counter;
    z *= 1.0f / counter;

    std::vector< float > dist(world.size(), 1e9);
    for(size_t i = 0; i < world.size(); ++i) {
        if(world[i].get_cell_type() == focal_cell_type) {
            if (!using_3d) dist[i] = (world[i].x_ - x) * (world[i].x_ - x) + (world[i].y_ - y) * (world[i].y_ - y);
            if (using_3d)  dist[i] = (world[i].x_ - x) * (world[i].x_ - x) + (world[i].y_ - y) * (world[i].y_ - y) + (world[i].z_ - z) * (world[i].z_ - z);
          }
     }

    auto min = std::min_element(dist.begin(), dist.end());
    return(std::distance(dist.begin(), min));
  }

  size_t find_central_cell(const std::vector< size_t >& positions) const {
    // first calculate average x and y of cell type:
    float x = 0.f;
    float y = 0.f;
    float z = 0.f;
    int counter = 0;
    for(const auto& i : positions) {
        x += world[i].x_;
        y += world[i].y_;
        if (using_3d) z += world[i].z_;
        counter++;
      }

    x *= 1.0f / counter;
    y *= 1.0f / counter;
    z *= 1.0f / counter;

    std::vector< float > dist(positions.size(), 1e9);
    for (size_t i = 0; i < positions.size(); ++i) {
        size_t pos = positions[i];
        if (!using_3d) dist[i] = (world[pos].x_ - x) * (world[pos].x_ - x) + (world[pos].y_ - y) * (world[pos].y_ - y);
        if (using_3d)  dist[i] = (world[pos].x_ - x) * (world[pos].x_ - x) + (world[pos].y_ - y) * (world[pos].y_ - y) + (world[pos].z_ - z) * (world[pos].z_ - z);
      }

    auto min = std::min_element(dist.begin(), dist.end());
    return(std::distance(dist.begin(), min));
  }

  void initialize_full() {
    for(auto& i : world) {
        change_cell_type(i.pos, normal);
      }

    parameters.initial_number_cancer_cells = static_cast<int>(0.1f * world.size());
    add_cells(cancer);

    // just for safety, do a full scan:
    for(auto& i : world) {
        update_growth_prob(i.pos);
        update_death_prob(i.pos);
      }

    infect_center(0.1f);

    // and update again, just for added safety:
    for(auto& i : world) {
        update_growth_prob(i.pos);
        update_death_prob(i.pos);
      }
  }

  void test_change_cell_type(const size_t& pos,
                             const cell_type& new_cell_type) {
    change_cell_type(pos, new_cell_type);
  }

  void test_event(size_t event) {
    do_event(event);
  }

  void test_update_rates() {
    update_rates();
  }

  size_t test_pick_event(const std::array<float, 12>& v, float s) {
    return pick_event(v, s);
  }

  void test_ask_infect_neighbours(size_t depth, size_t pos) {
    ask_infect_neighbours(depth, pos, 1);
  }
  void test_infect_periphery(float frac) {
    infect_periphery(frac);
  }
  void test_infect_random(float frac) {
    infect_random(frac);
  }
  void test_infect_center(float frac) {
    infect_center(frac);
  }

  void test_infect_center_largest(float frac) {
    infect_center_largest(frac);
  }
  void test_infect_all_cancer() {
    infect_all_cancer();
  }
  void test_infect_long_distance(size_t pos) {
    infect_long_distance(pos);
  }

  float get_rates(size_t event) {
    return rates[event];
  }

  void update_grow_params() override {
    parameters.time_adding_virus += parameters.time_adding_cancer;
    parameters.time_adding_virus_2 += parameters.time_adding_virus;
  }

  void change_cell_type(const size_t& pos,
                        const cell_type& new_cell_type) {
    cell_type previous_type = world[pos].get_cell_type();
    world[pos].set_cell_type(new_cell_type);

    update_death_prob(pos, previous_type, new_cell_type);

    update_growth_prob(pos);
    for(const auto& i : world[pos].neighbors) {
        update_growth_prob(i->pos);
      }
    update_count(previous_type, new_cell_type);
  }

  void initialize_network(std::vector< std::vector< voronoi_point > >& all_polys,
                          grid_type used_grid_type) override {
    // initialize default.
    for(size_t i = 0; i < 4; ++i) {
        growth_prob[i] = binned_distribution(sq_size, num_cells);
        death_prob[i] = binned_distribution(sq_size, num_cells);
    }

    for(auto& i : world) {
        i.prob_normal_infected = parameters.prob_normal_infection;
      }

    if(parameters.use_voronoi_grid == false) {
        for(auto& i : world) {
            i.update_neighbors(world, sq_size);
            change_cell_type(i.pos, empty);
          }
      }
    if(parameters.use_voronoi_grid == true) {

        setup_voronoi(all_polys, used_grid_type,
                      num_cells, sq_size, rndgen);
        for(size_t i = 0; i < num_cells; ++i) {
            world[i].inv_num_neighbors = 1.f / world[i].neighbors.size();
            update_growth_prob(i);
            update_death_prob(i);
          }
      }

    if(parameters.start_setup == grow || parameters.start_setup == converge) {
        add_cells(normal);
        for(size_t i = 0; i < num_cells; ++i) {
            update_growth_prob(i);
            update_death_prob(i);
          }
      }

    if(parameters.start_setup == full) {
        initialize_full();
    }

    for(size_t i = 0; i < growth_prob.size(); ++i) {
        growth_prob[i].update_all();
        death_prob[i].update_all();
    }
  }

private:


  std::array< float, 12> rates;

  std::vector<double> long_distance_infection_probability;

  void update_rates() {
    rates[0] = parameters.birth_normal   * growth_prob[normal].get_total_sum();
    rates[1] = parameters.death_normal   * death_prob[normal].get_total_sum();

    rates[2] = parameters.birth_cancer   * growth_prob[cancer].get_total_sum();
    rates[3] = parameters.death_cancer   * death_prob[cancer].get_total_sum();

    rates[4] = parameters.birth_infected * growth_prob[infected].get_total_sum();
    rates[5] = parameters.death_infected * death_prob[infected].get_total_sum();

    rates[6] = parameters.birth_cancer_resistant * growth_prob[resistant].get_total_sum();
    rates[7] = parameters.death_cancer_resistant * death_prob[resistant].get_total_sum();

    // t-cell related rates.
    rates[8]  = 0.f;
    rates[9]  = 0.f;
    rates[10] = 0.f;
    rates[11] = 0.f;

    if (parameters.t_cell_sensitivity[normal])       rates[8]  = parameters.t_cell_rate * t_cell_death_prob[normal].get_total_sum();
    if (parameters.t_cell_sensitivity[cancer])       rates[9]  = parameters.t_cell_rate * t_cell_death_prob[cancer].get_total_sum();
    if (parameters.t_cell_sensitivity[infected])     rates[10] = parameters.t_cell_rate * t_cell_death_prob[infected].get_total_sum() * parameters.t_cell_infected_relative_rate;
    if (parameters.t_cell_sensitivity[resistant])    rates[11] = parameters.t_cell_rate * t_cell_death_prob[resistant].get_total_sum();
  }

  size_t pick_event(const std::array< float, 12>& rates, float sum) {
    float r = rndgen.uniform() * sum;
    for(size_t i = 0; i < rates.size(); ++i) {
        r -= rates[i];
        if(r <= 0) {
            return i;
          }
      }
    return 0;
  }

  void do_event(size_t event) {

    switch(event) {
      case 0: {
          implement_growth(normal);       // birth normal
          break;
        }
      case 1: {
          implement_death(normal, false);        // death normal
          break;
        }

      case 2: {
          implement_growth(cancer);       // birth cancer
          break;
        }
      case 3: {
          implement_death(cancer, false);        // death cancer
          break;
        }

      case 4: {
          implement_growth(infected);     // birth infection
          break;
        }
      case 5: {
          implement_death(infected, false);     // death infection
          break;
        }

      case 6: {
          implement_growth(resistant);     // birth infection
          break;
        }
      case 7: {
          implement_death(resistant, false);     // death infection
          break;
        }

      case 8: {
          implement_death(normal, true); // normal death by t-cell
          break;
      }
      case 9: {
        implement_death(cancer, true); // cancer death by t-cell
        break;
      }
      case 10: {
        implement_death(infected, true); // infected death by t-cell
         break;
      }
      case 11: {
        implement_death(resistant, true); // resistant death by t-cell
        break;
      }

      default: {
          // do nothing
          break;
        }
      }
  }

  void implement_death(const cell_type& parent, bool by_t_cell) {
      size_t position_of_dying_cell;
      if (by_t_cell) {
          position_of_dying_cell = t_cell_death_prob[parent].draw_explicit(rndgen);
      } else {
          position_of_dying_cell = death_prob[parent].draw_explicit(rndgen);
      }
      change_cell_type(position_of_dying_cell, empty);

      if (parent == infected && parameters.prob_infection_upon_death > 0.f) {
          infect_long_distance(position_of_dying_cell);
      }
      if (parent == infected && parameters.t_cell_increase > 0) {
          increase_t_cell_concentration(position_of_dying_cell);
      }
  }

  void implement_growth(const cell_type& parent) {

    // drawing position of growth is slow/bottleneck:
    size_t position_of_grown_cell = growth_prob[parent].draw_explicit(rndgen);

    cell_type new_type = parent;

    if(parent == cancer) {
        if(rndgen.uniform() < parameters.freq_resistant) {
            new_type = resistant;
          }
      }

    change_cell_type(position_of_grown_cell, new_type);
  }


  void remove_entries(std::vector< size_t >& source,
                      const std::vector< size_t >& to_remove) {
    for(auto i : to_remove) {
        for(size_t j = 0; j < source.size(); ++j)  {
            if(source[j] == i) {
                source[j] = source.back();
                source.pop_back();
                break;
              }
          }
      }
  }

  void infect_random(float fraction) {
    size_t num_cancer_cells = num_cell_types[cancer];

    size_t infected_cells = 0;
    size_t to_be_infected = static_cast<size_t>(fraction * num_cancer_cells);
    if(to_be_infected == 0) return;

    std::vector< size_t > cancer_pos(num_cancer_cells);
    int j = 0;
    for(const auto& i : world) {
        if(i.get_cell_type() == cancer) {
            cancer_pos[j] = static_cast<size_t>(i.pos);
            j++;
          }
      }

    while(infected_cells < to_be_infected && num_cancer_cells > 0) {
        //size_t position_of_grown_cell = static_cast<size_t>(death_prob_rnd[cancer](rndgen.rndgen_));
        size_t index = static_cast<size_t>(rndgen.random_number(cancer_pos.size()));
        size_t position_of_infected_cell = cancer_pos[ index ];

        change_cell_type(position_of_infected_cell, infected);

        num_cancer_cells--; // easy count for now.
        infected_cells++;
        cancer_pos[index] = cancer_pos.back();
        cancer_pos.pop_back();
      }

    for(size_t i = 0; i < num_cells; ++i) {
        update_growth_prob(i);
        update_death_prob(i);
      }
  }

  void infect_center(float fraction) {
    // infect_center_largest();
    size_t num_cancer_cells = num_cell_types[cancer];

    size_t to_be_infected = static_cast<size_t>(fraction * num_cancer_cells);
    if(to_be_infected == 0) return;

    //now find starting cell to infect
    size_t focal_pos = find_central_cell(cancer);
    std::vector<size_t> cells_turned(1, focal_pos);

    change_cell_type(focal_pos, infected);

    size_t counter = 0;
    while(cells_turned.size() < to_be_infected && num_cancer_cells > 0) {
        focal_pos = cells_turned[counter];
        for(size_t i = 0; i < world[focal_pos].neighbors.size(); ++i) {
            size_t other_pos = world[focal_pos].neighbors[i]->pos;
            if(world[other_pos].get_cell_type() == cancer) {

                change_cell_type(other_pos, infected);
                cells_turned.push_back(other_pos);
                num_cancer_cells--;

                if(cells_turned.size() >= to_be_infected) break;
              }
          }
        counter++;
        if(counter > cells_turned.size()) break;
      }

    for(const auto& i : cells_turned) {
        update_growth_prob(i);
        update_death_prob(i);
      }
  }


  void infect_center_largest(float fraction) {
    size_t num_cancer_cells = num_cell_types[cancer];
    if (num_cancer_cells == 0) return;
    // here, we inject cancer in the largest tumour mass.
    // thus, we first have to collect all tumour cells, and group them.

    std::vector< size_t > cancer_cell_pos(num_cancer_cells);
    size_t cnt = 0;
    size_t j = 0;
    for (const auto& i : world) {
        if (i.get_cell_type() == cancer) {
            cancer_cell_pos[cnt] = j;
            cnt++;
          }
        j++;
      }

    std::vector< std::vector< size_t > > clusters;
    std::vector< size_t > cluster_sizes;
    while(!cancer_cell_pos.empty()) {
        std::vector< size_t > cluster;
        cluster.push_back(cancer_cell_pos[0]);

        for(size_t i = 0; i < cluster.size(); ++i) {
            std::vector< size_t > neighbours = world[cluster[i]].get_cancer_neighbours();
            for(size_t j = 0; j < neighbours.size(); ++j) {

                // this is not the fastest method, but it is very reliable.
                // probably a method using an unsorted set could work much faster.
                if(std::find(cluster.begin(), cluster.end(), neighbours[j]) ==
                   cluster.end()) {
                    cluster.push_back(neighbours[j]);
                  }
              }
          }
        // now we have to remove these entries
        remove_entries(cancer_cell_pos, cluster);
        cluster_sizes.push_back(cluster.size());
        clusters.push_back(cluster);
      }

    auto m = std::max_element(cluster_sizes.begin(), cluster_sizes.end());

    size_t largest_cluster = 0;
    for (size_t i = 0; i < cluster_sizes.size(); ++i) {
        if (cluster_sizes[i] == *m) {
            largest_cluster = i;
            break;
          }
      }
    size_t to_be_infected = static_cast<size_t>(fraction * *m);
    if(to_be_infected == 0) return;

    if(*m < to_be_infected) {
        to_be_infected = *m;
      }

    size_t index_central_cell = find_central_cell(clusters[largest_cluster]);

    size_t focal_pos = clusters[largest_cluster][index_central_cell];
    std::vector<size_t> cells_turned(1, focal_pos);

    change_cell_type(focal_pos, infected);

    size_t counter = 0;
    while(cells_turned.size() < to_be_infected && num_cancer_cells > 0) {
        focal_pos = cells_turned[counter];
        for(size_t i = 0; i < world[focal_pos].neighbors.size(); ++i) {
            size_t other_pos = world[focal_pos].neighbors[i]->pos;
            if(world[other_pos].get_cell_type() == cancer) {

                change_cell_type(other_pos, infected);
                cells_turned.push_back(other_pos);
                num_cancer_cells--;

                if(cells_turned.size() >= to_be_infected) break;
              }
          }
        counter++;
        if(counter > cells_turned.size()) break;
      }

    for(const auto& i : cells_turned) {
        update_growth_prob(i);
        update_death_prob(i);
      }
  }

  void infect_periphery(float fraction) {
    size_t num_cancer_cells = num_cell_types[cancer];
    if (num_cancer_cells == 0) return;
    // here, we inject cancer in the largest tumour mass.
    // thus, we first have to collect all tumour cells, and group them.

    std::vector< size_t > cancer_cell_pos(num_cancer_cells);
    size_t cnt = 0;
    size_t j = 0;
    for (const auto& i : world) {
        if (i.get_cell_type() == cancer) {
            cancer_cell_pos[cnt] = j;
            cnt++;
          }
        j++;
      }

    std::vector< std::vector< size_t > > clusters;
    std::vector< size_t > cluster_sizes;
    while(!cancer_cell_pos.empty()) {
        std::vector< size_t > cluster;
        cluster.push_back(cancer_cell_pos[0]);

        for(size_t i = 0; i < cluster.size(); ++i) {
            std::vector< size_t > neighbours = world[cluster[i]].get_cancer_neighbours();
            for(size_t j = 0; j < neighbours.size(); ++j) {

                // this is not the fastest method, but it is very reliable.
                // probably a method using an unsorted set could work much faster.
                if(std::find(cluster.begin(), cluster.end(), neighbours[j]) ==
                   cluster.end()) {
                    cluster.push_back(neighbours[j]);
                  }
              }
          }
        // now we have to remove these entries
        remove_entries(cancer_cell_pos, cluster);
        cluster_sizes.push_back(cluster.size());
        clusters.push_back(cluster);
      }

    auto m = std::max_element(cluster_sizes.begin(), cluster_sizes.end());

    size_t largest_cluster = 0;
    for (size_t i = 0; i < cluster_sizes.size(); ++i) {
        if (cluster_sizes[i] == *m) {
            largest_cluster = i;
            break;
          }
      }

    // now we have the largest cluster.
    // let's now determine the distance to the center of the cluster.

    // first, we determine the center:
    float x = 0.f;
    float y = 0.f;
    int counter = 0;
    for(const auto& i : clusters[largest_cluster]) {
        x += world[i].x_;
        y += world[i].y_;
        counter++;
      }

    x *= 1.0f / counter;
    y *= 1.0f / counter;

    struct peri_cell {
      size_t pos;
      float freq_cancer;
      float dist_to_center;
    };

    std::vector< peri_cell > periphery;

    for (size_t i = 0; i < clusters[largest_cluster].size(); ++i) {
        size_t pos = clusters[largest_cluster][i];

        float freq_cancer = world[pos].freq_type_neighbours(cancer);

        if(freq_cancer < 1.0f &&
           world[pos].get_cell_type() == cancer) {
            peri_cell add;
            add.dist_to_center = (world[pos].x_ - x) * (world[pos].x_ - x) +
                (world[pos].y_ - y) * (world[pos].y_ - y);
            add.pos = pos;
            add.freq_cancer = freq_cancer;
            periphery.push_back(add);
          }
      }

    size_t to_be_infected = static_cast<size_t>(fraction * periphery.size());
    if(to_be_infected == 0) return;

    std::sort(periphery.begin(), periphery.end(),
              [](auto const& a, auto const& b) {
        if (a.dist_to_center == b.dist_to_center) {
            return a.freq_cancer < b.freq_cancer;
          }
        return a.dist_to_center > b.dist_to_center;});


    std::vector< size_t > cells_turned;
    for (size_t i = 0; i < to_be_infected; ++i) {
        size_t focal_pos = periphery[i].pos;
        if (world[focal_pos].get_cell_type() == cancer) {
            change_cell_type(focal_pos, infected);
            cells_turned.push_back(focal_pos);
          }
      }

    for(const auto& i : cells_turned) {
        update_growth_prob(i);
        update_death_prob(i);
      }
  }
  void infect_all_cancer() {
    for(auto& i : world) {
        if(i.get_cell_type() == cancer) {
            change_cell_type(i.pos, infected);
          }
      }

    for(const auto& i : world) {
        update_growth_prob(i.pos);
        update_death_prob(i.pos);
      }
  }

  void infect_long_distance(size_t pos) {
    size_t identifier = rndgen.random_number(1e10);
    // identifier is used as a unique ID.
    ask_infect_neighbours(parameters.distance_infection_upon_death,
                          pos, identifier);
  }

  void ask_infect_neighbours(size_t depth,
                             size_t pos,
                             size_t identifier) {

    float p = long_distance_infection_probability[depth];

    if (depth == 0) return;

    depth--;
    for(auto& n : world[pos].neighbors) {

        if (n->get_identifier() != identifier) {
            if(n->get_cell_type() == cancer) {
                if(rndgen.uniform() < p) {
                    change_cell_type(n->pos, infected);
                  }
              }
            if(n->get_cell_type() == normal) {
                if(rndgen.uniform() < p) {
                    if(rndgen.uniform() < n->prob_normal_infected) {
                        change_cell_type(n->pos, infected);
                      }
                  }
              }
            n->set_identifier(identifier);
          }
        ask_infect_neighbours(depth, n->pos, identifier);
      }
  }

  void update_count(cell_type old_type, cell_type new_type) {
    num_cell_types[old_type]--;
    num_cell_types[new_type]++;
  }

  void setup_voronoi(std::vector< std::vector< voronoi_point > >& all_polys,
                     grid_type used_grid_type,
                     size_t num_cells,
                     size_t sq_size,
                     rnd_t& rndgen) {

    // std::cout << "Generating centre points\n";
    std::vector< voronoi_point > v(num_cells);


    // we make regular grid
    for(size_t i = 0; i < num_cells; ++i) {

        float x, y;

        if (used_grid_type == grid_type::hexagonal) {
            x = i % sq_size;
            y = i / sq_size;
            if ((i / sq_size) % 2 == 0) x += 0.5f;
          } else {
            x = rndgen.uniform() * sq_size;
            y = rndgen.uniform() * sq_size;
          }

        v[i] = voronoi_point(x, y);
      }

    cinekine::voronoi::Sites sites;
    // << "converting centre points to vertices\n";
    for(auto i : v) {
        cinekine::voronoi::Vertex temp_vertex(i.x_, i.y_);
        sites.push_back(temp_vertex);
      }

    // std::cout << "creating voronoi graph\n";
    cinekine::voronoi::Graph graph = cinekine::voronoi::build(std::move(sites), sq_size, sq_size);

    std::vector< std::vector< voronoi_edge > > all_edges(world.size());

    // std::cout << "Ready to build world\n";
    // std::cout << "collecting all edges\n";
    for(const auto& cell : graph.cells()) {

        size_t site_index = static_cast<size_t>(cell.site);

        cinekine::voronoi::Site focal_site = graph.sites()[site_index];

        world[site_index].x_ = focal_site.x;
        world[site_index].y_ = focal_site.y;

        for(const auto& edge : cell.halfEdges) {
            cinekine::voronoi::Edge focal_edge = graph.edges()[edge.edge];
            voronoi_point start(focal_edge.p0.x, focal_edge.p0.y);
            voronoi_point end(  focal_edge.p1.x, focal_edge.p1.y);

            voronoi_edge local_edge(start, end, focal_edge.leftSite, focal_edge.rightSite);

            if(local_edge.calc_dist() > 1e-2) {
                all_edges[site_index].push_back(local_edge);
              }
          }
      }

    // std::cout << "implementing all edges\n";
    for(auto i : all_edges) {
        for(auto edge : i) {
            size_t left  = edge.left;
            size_t right = edge.right;

            if(left < world.size() && right < world.size()) {
                world[left].template add_neighbor<NODE>(world, right);
                world[right].template add_neighbor<NODE>(world, left);
              }
          }
      }

    // std::cout << "clean edges and add polies for plotting\n";
    for(size_t i = 0; i < num_cells; ++i) {
        std::vector< voronoi_point > poly = voronoi_tools::clean_edges(all_edges[i], i);
        all_polys.push_back(poly);
      }
  }


};

inline std::unique_ptr<simulation> create_simulation(bool use_3d,
                                                     const Param& p) {
  if (!use_3d) {
    return std::unique_ptr<simulation>(new simulation_impl<node_2d>(p));
  }
  else {
    return std::unique_ptr<simulation>(new simulation_impl<node_3d>(p, true));
  }
}

inline std::unique_ptr<simulation> create_simulation(const Param& p) {
  if (!p.using_3d) {
    return std::unique_ptr<simulation>(new simulation_impl<node_2d>(p));
  }
  else {
    return std::unique_ptr<simulation>(new simulation_impl<node_3d>(p, true));
  }
}




#endif /* simulation_hpp */
