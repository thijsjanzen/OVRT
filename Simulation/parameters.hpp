//
//  parameters.h
//
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#ifndef parameters_h
#define parameters_h

#include <cstdio>
#include <cmath>
#include <random>

enum infection_routine {random_infection, center_infection, periphery_infection};
enum start_type {full, grow, converge, empty_grid};
enum grid_type {regular, voronoi, hexagonal};
enum t_cell_type {cancer_cell, infected_cell, both};


struct Param {

  size_t seed;

  int maximum_time;

  float birth_normal;
  float birth_cancer;
  float birth_infected;
  float death_normal;
  float death_cancer;
  float death_infected;

  float birth_cancer_resistant;
  float death_cancer_resistant;
  float freq_resistant;


  int time_adding_virus;
  int time_adding_virus_2;
  int time_adding_cancer;

  size_t initial_number_cancer_cells;
  size_t initial_number_normal_cells;

  float percent_infected;
  float percent_infected_2;
  infection_routine infection_type;
  infection_routine infection_type_2;

  float prob_normal_infection;

  float distance_infection_upon_death;
  float prob_infection_upon_death;

  size_t sq_num_cells;
  size_t sq_num_pixels;
  bool use_voronoi_grid;

  bool using_3d;

  start_type start_setup;
  t_cell_type t_cell_effect;

  float t_cell_increase;
  float evaporation;
  float diffusion;
  float t_cell_rate;
  float t_cell_threshold;
  float t_cell_density_scaler;
  float t_cell_inflection_point;

  Param() {
    std::random_device rd;
    seed = rd();

    // default values, these don't necessarily make sense!
    maximum_time = 1000; // default from Berg et al 2019
    time_adding_virus = 100; // default from Berg et al 2019
    time_adding_cancer = 500; // normally the population is first allowed to stabilize

    infection_type = random_infection;
    percent_infected = 0.1f;

    time_adding_virus_2 = 500000; // extremely high value to make it unlikely to trigger by default
    infection_type_2 = random_infection;
    percent_infected_2 = 0.1f;

    initial_number_cancer_cells = 500;
    initial_number_normal_cells = 1000;

    birth_normal = 0.5f;
    birth_cancer = 1.0f;
    birth_infected = 1.2f;
    death_normal = 0.2f;
    death_cancer = 0.1f;
    death_infected = 0.1f;

    death_cancer_resistant = 0.1f;
    birth_cancer_resistant = 0.9f;

    prob_normal_infection = 0.0f;
    freq_resistant = 0.0f;
    distance_infection_upon_death = 1;
    prob_infection_upon_death = 0.f;

    start_setup = converge;

    sq_num_cells = 100;
    sq_num_pixels = 100;

    use_voronoi_grid = false;
    using_3d = false;


    t_cell_effect = both;
    t_cell_increase = 0.f;
    evaporation = 0.01f;
    diffusion = 0.1f;
    t_cell_rate = 10;
    t_cell_threshold = 0.2f;
    t_cell_density_scaler = 1.0f;
  }
};

#endif /* parameters_h */
