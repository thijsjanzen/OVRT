#include <cstring>
#include <chrono>

#include "../Simulation/simulation.hpp"
#include "../Simulation/analysis.hpp"
#include "config_parser.h"

// forward declaration
bool file_exists (const std::string& name);
void read_parameters_from_ini(Param& p, const std::string file_name);
void obtain_equilibrium(simulation& Simulation, const Param& all_parameters);

int main(int argc, char *argv[]) {

  std::cout << "Welcome to this In Silico Simulation of oncolytic tumor virotherapy\n";
  std::cout << "Copyright 2019 - 2024, D. Bhatt, T. Janzen & F.J. Weissing\n";
  std::cout << "This is version: 0.9.6\n";

  std::cout << "All files are to be found in this folder: \n";
  std::cout << argv[0] << "\n";

  InputParser input(argc, argv);

  std::string file_name = "config.ini";

  const std::string &filename = input.getCmdOption("-f");
  if (!filename.empty()){
          file_name = filename;
  }

  Param all_parameters;

  read_parameters_from_ini(all_parameters, file_name);

  read_from_command_line(input, all_parameters);

  float max_t = 0;
  std::array<size_t, 5> cell_counts = do_analysis(all_parameters, max_t);
  std::string outcome = get_outcome(cell_counts);

  std::cout << outcome << "\n";

  if(!file_exists("output.txt")) {
      // write header to file
      std::ofstream outfile("output.txt");
      outfile << "birth_virus"         << "\t"
              << "death_virus"         << "\t"
              << "birth_cancer"        << "\t"
              << "death_cancer"        << "\t"
              << "freq_resistant"      << "\t"
              << "outcome"             << "\t"
              << "num_normal_cells"    << "\t"
              << "num_cancer_cells"    << "\t"
              << "num_infected_cells"  << "\t"
              << "num_resistant_cells" << "\t"
              << "num_empty_cells"     << "\t"
              << "total_runtime"       << "\n";
      outfile.close();
   }

  std::ofstream outfile("output.txt", std::ios::app);
  outfile << all_parameters.birth_infected << "\t"
          << all_parameters.death_infected << "\t"
          << all_parameters.birth_cancer   << "\t"
          << all_parameters.death_cancer   << "\t"
          << all_parameters.freq_resistant << "\t"
          << outcome                       << "\t";
  for(size_t i = 0; i < 5; ++i) {
      outfile << cell_counts[i] << "\t";
  }
  outfile << max_t << "\t";
  outfile << "\n";
  outfile.close();

  return 0;
}

void read_parameters_from_ini(Param& p, const std::string file_name) {

  ConfigFile from_config(file_name);

  p.seed = from_config.getValueOfKey<size_t>("seed");

  p.maximum_time = from_config.getValueOfKey<int>("maximum_time");
  p.time_adding_cancer = from_config.getValueOfKey<int>("time_adding_cancer");
  p.time_adding_virus  = from_config.getValueOfKey<int>("time_adding_virus");

  p.initial_number_cancer_cells = from_config.getValueOfKey<size_t>("initial_number_cancer_cells");
  p.initial_number_normal_cells = from_config.getValueOfKey<size_t>("initial_number_normal_cells");

  p.birth_normal = from_config.getValueOfKey<float>("birth_normal");
  p.death_normal = from_config.getValueOfKey<float>("death_normal");

  p.birth_cancer = from_config.getValueOfKey<float>("birth_cancer");
  p.death_cancer = from_config.getValueOfKey<float>("death_cancer");

  p.birth_infected = from_config.getValueOfKey<float>("birth_infected");
  p.death_infected = from_config.getValueOfKey<float>("death_infected");

  p.birth_cancer_resistant = from_config.getValueOfKey<float>("birth_resistant");
  p.death_cancer_resistant = from_config.getValueOfKey<float>("death_resistant");

  p.percent_infected = from_config.getValueOfKey<float>("percent_infected");
  p.prob_normal_infection = from_config.getValueOfKey<float>("prob_normal_infection");
  p.freq_resistant = from_config.getValueOfKey<float>("freq_resistant");

  p.distance_infection_upon_death = from_config.getValueOfKey<float>("distance_infection_upon_death");
  p.prob_infection_upon_death = from_config.getValueOfKey<float>("prob_infection_upon_death");

  p.sq_num_cells = from_config.getValueOfKey<size_t>("sq_num_cells");

  p.using_3d = from_config.getValueOfKey<bool>("using_3d");

  p.infection_type = random_infection;
  auto infection_string = from_config.getValueOfKey<std::string>("infection_type");
  if(infection_string == "Random")
    p.infection_type = random_infection;
  if(infection_string == "Center")
    p.infection_type = center_infection;

  auto start_string = from_config.getValueOfKey<std::string>("start_setup");
  if(start_string == "Grow")
    p.start_setup = grow;
  if(start_string == "Full")
    p.start_setup = full;
  if(start_string == "Random")
    p.start_setup = random_grid;

  auto grid_string = from_config.getValueOfKey<std::string>("grid_type");
  if(grid_string == "regular")
    p.use_voronoi_grid = false;
  if(grid_string == "voronoi")
    p.use_voronoi_grid = true;

  auto t_cell_sensitivity_stromal_string = from_config.getValueOfKey<std::string>("t_cell_stromal");
  if (t_cell_sensitivity_stromal_string == "sensitive") {
      p.t_cell_sensitivity[normal] = true;
  } else {
      p.t_cell_sensitivity[normal] = false;
  }

  auto t_cell_sensitivity_cancer_string = from_config.getValueOfKey<std::string>("t_cell_cancer");
  if (t_cell_sensitivity_cancer_string == "sensitive") {
      p.t_cell_sensitivity[cancer] = true;
  } else {
      p.t_cell_sensitivity[cancer] = false;
  }

  auto t_cell_sensitivity_infected_string = from_config.getValueOfKey<std::string>("t_cell_infected");
  if (t_cell_sensitivity_infected_string == "sensitive") {
      p.t_cell_sensitivity[infected] = true;
  } else {
      p.t_cell_sensitivity[infected] = false;
  }

  auto t_cell_sensitivity_resistant_string = from_config.getValueOfKey<std::string>("t_cell_resistant");
  if (t_cell_sensitivity_resistant_string == "sensitive") {
      p.t_cell_sensitivity[resistant] = true;
  } else {
      p.t_cell_sensitivity[resistant] = false;
  }

  p.t_cell_infected_relative_rate = from_config.getValueOfKey<float>("relative_rate_infected");

  return;
}

bool file_exists (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}
