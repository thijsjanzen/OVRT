#include "analysis.hpp"
#include <fstream>

std::array<size_t, 5> do_analysis(Param all_parameters) {

  std::unique_ptr<simulation> Simulation = create_simulation(all_parameters.using_3d, all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  if (all_parameters.use_voronoi_grid == false) {
    Simulation->initialize_network(filler, regular);
  } else {
      Simulation->initialize_network(filler, voronoi);
  }


  Simulation->t = 0.f;
  float prev_t = Simulation->t;
  std::array<size_t, 5> cell_counts;

  if(all_parameters.start_setup == converge) {
      std::cout << "simulating until having reached equilibrium with normal cells\n";
      Simulation->obtain_equilibrium(false); // this obtains a fully grown grid, with normal cells
  }

  Simulation->t = 0.f;
  bool cancer_added = false;
  bool virus_added = false;
  bool virus_added_2 = false;

  auto prev_timepoint = std::chrono::steady_clock::now();
  auto start_t = prev_timepoint;

  float prev_cast_t = static_cast<int>(Simulation->t);

  while(Simulation->t < all_parameters.maximum_time) {
          Simulation->update_one_step();

          if(prev_t < all_parameters.time_adding_cancer &&
             Simulation->t >= all_parameters.time_adding_cancer &&
             cancer_added == false) {
           //   std::cout << "adding cancer!\n";
              Simulation->add_cells(cancer);
              Simulation->t = 0.f; // reset time
              prev_cast_t = static_cast<int>(Simulation->t);
              cancer_added = true;
          }
          if(prev_t < all_parameters.time_adding_virus &&
             Simulation->t >= all_parameters.time_adding_virus &&
             virus_added   == false  &&
             cancer_added  == true) {
           //  std::cout << "adding virus!\n";
             Simulation->add_infected(all_parameters.infection_type,
                                     all_parameters.percent_infected);
             Simulation->t = 0.f;
             prev_cast_t = static_cast<int>(Simulation->t);
             virus_added = true;
          }

          if(prev_t < all_parameters.time_adding_virus_2 &&
             Simulation->t >= all_parameters.time_adding_virus_2 &&
             virus_added == true &&
             virus_added_2 == false &&
             cancer_added  == true) {
           //   std::cout << "adding virus for the second time!\n";
             Simulation->add_infected(all_parameters.infection_type_2,
                                     all_parameters.percent_infected_2);
             Simulation->t = 0.f;
             prev_cast_t = static_cast<int>(Simulation->t);
             virus_added_2 = true;
          }

          auto next_t = std::chrono::steady_clock::now();

          // auto diff_t = std::chrono::duration_cast<std::chrono::seconds>(next_t - prev_timepoint).count();

          int cast_t = static_cast<int>(Simulation->t);

          if((cast_t - prev_cast_t) >= 10) {
              prev_cast_t = cast_t;

            cell_counts = Simulation->get_count_cell_types();  //Simulation->num_cell_types;

              if(all_parameters.start_setup == grow ||
                 all_parameters.start_setup == converge) {

                if(cell_counts[cancer] < 1 && virus_added == true) {
                  //  std::cout << "cancer eradicated\n";
                    break;
                }
                if(cell_counts[normal] < 1 && virus_added == true) {
                  //  std::cout << "normal tissue gone\n";
                    break;
                }
                if(cell_counts[infected] < 1 && virus_added == true &&
                   cell_counts[normal] < 1) {
                  //  std::cout << "virus wiped out\n";
                    break;
                }
              }
              if(all_parameters.start_setup == full) {
                  if(cell_counts[cancer] < 1) {
                   //   std::cout << "cancer eradicated";
                      break;
                  }
                  if(cell_counts[normal] < 1) {
                    //  std::cout << "normal tissue gone\n";
                      break;
                  }
               }

              prev_timepoint = next_t;
            }

          prev_t = Simulation->t;
  }

  // we want to know:
  // A: tumor eradication, e.g. only normal cells remain
  // B: tumor victory, e.g. only tumor cells remain
  // C: co-existence of the three populations
  cell_counts = Simulation->get_count_cell_types();
  return cell_counts;
}

std::string get_outcome(const std::array<size_t, 5>& cell_counts) {
  std::array<float, 4> freq = {0.f, 0.f, 0.f, 0.f};
  freq[normal] = cell_counts[normal];
  freq[cancer] = cell_counts[cancer];
  freq[infected] = cell_counts[infected];
  freq[resistant] = cell_counts[resistant];
  float total_num_cells = cell_counts[normal] + cell_counts[cancer] + cell_counts[infected] + cell_counts[resistant];
  // we want to know:
  // A: tumor eradication, e.g. only normal cells remain
  // B: tumor victory, e.g. only tumor cells remain
  // C: co-existence of the three populations
  // D: tumor & resistant cells remain, no virus or normal cells.
  for(size_t i = 0; i < 4; ++i) freq[i] *= 1.0f / total_num_cells;

  if(freq[resistant] < 1e-6f) {

    if(freq[cancer] <= 1e-6f) {
        return "A";
    }

    if(freq[normal] <= 1e-6f     && freq[cancer] >= (1-1e-6f) && freq[infected] <= 1e-6f) {
        return "B";
    }

    return "C";
  } else {
    // frequency resistant is non-zero
    if(freq[infected] < freq[cancer]) {
        return "D";
    } else {
        return "C";
    }
  }

}
