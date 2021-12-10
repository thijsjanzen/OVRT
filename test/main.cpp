#define CATCH_CONFIG_MAIN
#include "catch.h"
#include "../Simulation/simulation.hpp"
#include "../Simulation/analysis.hpp"


TEST_CASE( "birth_death" )
{
  std::cout << "testing birth and death\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_2d> Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  Simulation.t = 0.f;

  // count cell types
  std::array<size_t, 5> cells = Simulation.count_cell_types();
  for(size_t i = 0; i < 4; ++i) {
     REQUIRE(cells[i] == 0);
  }
  REQUIRE(cells[4] == all_parameters.sq_num_cells *
                              all_parameters.sq_num_cells);


  // test update_rates

  Simulation.test_update_rates();
  for(auto c : {normal, cancer, infected, resistant}) {
   REQUIRE(Simulation.death_prob[c].get_total_sum() == 0.0);
   REQUIRE(Simulation.growth_prob[c].get_total_sum() == 0.0);
   REQUIRE(Simulation.get_rates(c) == 0.0);
  }

  Simulation.test_change_cell_type(0, normal);
  Simulation.test_update_rates();
  REQUIRE(Simulation.get_rates(1) == all_parameters.death_normal);

  Simulation.test_change_cell_type(3, cancer);
  Simulation.test_update_rates();
  REQUIRE(Simulation.get_rates(3) == all_parameters.death_cancer);

  Simulation.test_change_cell_type(0, empty);
  Simulation.test_change_cell_type(3, empty);
  Simulation.test_update_rates();




  // test change_chell_type
  // add a normal cell, then kill it

  for(auto c : {normal, cancer, infected, resistant}) {
    Simulation.test_change_cell_type(5000, c);

    REQUIRE(Simulation.world[5000].get_cell_type() ==
                      c);

    REQUIRE(Simulation.death_prob[c].get_value(5000) ==
                      1.0);

    Simulation.test_change_cell_type(5000, empty);

    REQUIRE(Simulation.world[5000].get_cell_type() ==
                      empty);

    REQUIRE(Simulation.death_prob[c].get_value(5000) ==
                      0.0);
   }

  // test do_event

  std::vector< cell_type > v = {normal, empty, cancer, empty,
                                infected, empty, resistant, empty};

   for(size_t i = 0; i < 8; ++i) {
    Simulation.test_event(i);

    REQUIRE(Simulation.world[0].get_cell_type() ==
                      v[i]);

    if (v[i] != empty) {
        REQUIRE(Simulation.death_prob[v[i]].get_value(0) ==
                          1.0);
    }
   }

   // test pick event
   std::array<float, 8 > vx = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
   REQUIRE(Simulation.test_pick_event(vx, 0.0) == 0);
   for(size_t i = 0; i < 8; ++i) {
    vx[i] = 1.0f;
    REQUIRE(Simulation.test_pick_event(vx, 1.0) == i);
    vx[i] = 0.0f;
   }
}

TEST_CASE( "check_outcome" )
{
  std::cout << "testing check_outcome\n";
  std::array<size_t, 5> cell_counts_normal = {100000, 0, 0, 0, 0};

  REQUIRE(get_outcome(cell_counts_normal) == "A");

  std::array<size_t, 5> cell_counts_cancer = {0, 1000000, 0, 0, 0};

  REQUIRE(get_outcome(cell_counts_cancer) == "B");

  std::array<size_t, 5> cell_counts_resistant = {0, 0, 0, 1000000, 0};

  REQUIRE(get_outcome(cell_counts_resistant) == "C");
}


TEST_CASE( "find_central_cell" )
{
  std::cout << "testing finding center cells\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_2d> Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  // first we try to find the central empty cell (50, 50), coordinate: 5000
  auto index = Simulation.find_central_cell(empty);
  auto x = Simulation.world[index].x_;
  auto y = Simulation.world[index].y_;
  REQUIRE(x == 49);
  REQUIRE(y == 49);

  size_t row_size = all_parameters.sq_num_cells;
  // now we add a square of cancer cells, and calculate the center
  std::vector< size_t > positions;
  for(size_t x = 10; x < 21; ++x) {
      for(size_t y = 10; y < 21; ++y) {
          size_t pos = y * row_size + x;
          Simulation.test_change_cell_type(pos, cancer);
          positions.push_back(pos);
        }
  }
  auto index2 = Simulation.find_central_cell(cancer);
  auto x2 = Simulation.world[index2].x_;
  auto y2 = Simulation.world[index2].y_;
  REQUIRE(x2 == 15);
  REQUIRE(y2 == 15);

  auto index3 = Simulation.find_central_cell(positions);
  auto x3 = Simulation.world[positions[index3]].x_;
  auto y3 = Simulation.world[positions[index3]].y_;
  REQUIRE(x2 == x3);
  REQUIRE(y2 == y3);
}

TEST_CASE( "add_cells" )
{
  std::cout << "testing adding cells\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;
  all_parameters.initial_number_normal_cells = 1000;
  all_parameters.initial_number_cancer_cells = 100;

  simulation_impl<node_2d> Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  std::array<size_t, 5> cell_cnt = Simulation.count_cell_types();
  REQUIRE(cell_cnt[normal] == 0);

  Simulation.add_cells(normal);
  cell_cnt = Simulation.count_cell_types();
  REQUIRE(cell_cnt[normal] ==
                    all_parameters.initial_number_normal_cells);

  Simulation.add_cells(cancer);
  cell_cnt = Simulation.count_cell_types();
  REQUIRE(cell_cnt[cancer] ==
                    all_parameters.initial_number_cancer_cells);
}

TEST_CASE( "setup_types")
{
  std::cout << "testing setup types\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = full;

  simulation_impl<node_2d> Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);
  std::array<size_t, 5> cell_cnt = Simulation.count_cell_types();

  size_t total_num_cells = all_parameters.sq_num_cells *
                           all_parameters.sq_num_cells;

  REQUIRE(cell_cnt[normal] ==
                    total_num_cells * 0.9);
  REQUIRE(cell_cnt[cancer] ==
                    total_num_cells * 0.09);
  REQUIRE(cell_cnt[infected] ==
                    total_num_cells * 0.01);
}

TEST_CASE( "ask_infect_neighbours")
{
  std::cout << "testing ask infect neighbours\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = full;
  all_parameters.distance_infection_upon_death = 1.0;
  all_parameters.prob_infection_upon_death = 1.f;

  simulation_impl<node_2d> Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  size_t row_size = all_parameters.sq_num_cells;
  size_t x = 50;
  for(size_t y = 40; y < 61; ++y) {
      size_t pos = y * row_size + x;
      Simulation.test_change_cell_type(pos, cancer);
  }
  size_t y2 = 50;
  for(size_t x2 = 40; x2 < 61; ++x2) {
      size_t pos = y2 * row_size + x2;
      Simulation.test_change_cell_type(pos, cancer);
  }

  // now we have a cross of individuals
  size_t initial_pos = 50 * 100 + 50;
  Simulation.test_change_cell_type(initial_pos, infected);
  Simulation.test_ask_infect_neighbours(all_parameters.distance_infection_upon_death,
                                        initial_pos);
  for(size_t x = 40; x < 61; ++x) {
    size_t pos = 50 * row_size + x;
    int dist_x = static_cast<int>(x) - static_cast<int>(50);
    if (dist_x < 0) dist_x *= -1;
    if (dist_x == 0) {
        continue;
    }

    auto ct = Simulation.world[pos].get_cell_type();

    if (dist_x <= static_cast<int>(all_parameters.distance_infection_upon_death)) {
      REQUIRE(ct == infected);
    } else {
      REQUIRE(ct == cancer);
    }

  }
}

TEST_CASE( "ask_infect_neighbours_two")
{
  std::cout << "testing ask infect neighbours\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = full;
  all_parameters.distance_infection_upon_death = 2.0;
  all_parameters.prob_infection_upon_death = 1.f;

  simulation_impl<node_2d> Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  size_t row_size = all_parameters.sq_num_cells;
  size_t x = 50;
  for(size_t y = 40; y < 61; ++y) {
      size_t pos = y * row_size + x;
      Simulation.test_change_cell_type(pos, cancer);
  }
  size_t y2 = 50;
  for(size_t x2 = 40; x2 < 61; ++x2) {
      size_t pos = y2 * row_size + x2;
      Simulation.test_change_cell_type(pos, cancer);
  }

  // now we have a cross of individuals
  size_t initial_pos = 50 * 100 + 50;
  Simulation.test_change_cell_type(initial_pos, infected);
  Simulation.test_ask_infect_neighbours(all_parameters.distance_infection_upon_death,
                                        initial_pos);
  for(size_t x = 40; x < 61; ++x) {
    size_t pos = 50 * row_size + x;
    int dist_x = static_cast<int>(x) - static_cast<int>(50);
    if (dist_x < 0) dist_x *= -1;
    if (dist_x == 0) {
        continue;
    }

    auto ct = Simulation.world[pos].get_cell_type();

    if (dist_x <= static_cast<int>(all_parameters.distance_infection_upon_death)) {
      REQUIRE(ct == infected);
    } else {
      REQUIRE(ct == cancer);
    }
  }

  for(size_t y = 40; y < 61; ++y) {
    size_t pos = 50 * row_size + y;
    int dist_y = static_cast<int>(y) - static_cast<int>(50);
    if (dist_y < 0) dist_y *= -1;
    if (dist_y == 0) {
        continue;
    }

    auto ct = Simulation.world[pos].get_cell_type();

    if (dist_y <= static_cast<int>(all_parameters.distance_infection_upon_death)) {
      REQUIRE(ct == infected);
    } else {
      REQUIRE(ct == cancer);
    }
  }


}

TEST_CASE( "random_stuff" )
{
  std::cout << "testing randomizer\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 10;

  simulation_impl<node_2d> Simulation(all_parameters);

  // test random numbers
  Simulation.growth_prob[0].update_entry(10, 1.0f);
  auto x = Simulation.growth_prob[0].draw_explicit(Simulation.rndgen);
  REQUIRE(x == 10);

  Simulation.growth_prob[0].update_entry(10, 0.0f);
  Simulation.growth_prob[0].update_entry(50, 0.01f);
  auto x2 = Simulation.growth_prob[0].draw_explicit(Simulation.rndgen);
  REQUIRE(x2 == 50);

  // now we populate with equal numbers
  int total_num_cells = all_parameters.sq_num_cells * all_parameters.sq_num_cells;
  for(int i = 0; i < total_num_cells; ++i) {
    Simulation.growth_prob[0].update_entry(i, 1.0f);
  }
  // and now we draw many numbers:
  int num_zero = 0;
  int num_ten  = 0;
  int num_repl = 100000;
  for(int r = 0; r < num_repl; ++r) {
      auto xx = Simulation.growth_prob[0].draw_explicit(Simulation.rndgen);
      if(xx == 0) num_zero++;
      if(xx == 10) num_ten++;
  }

  float expected_freq = 1.f / total_num_cells;
  float freq_zero     = 1.f * num_zero / num_repl;
  float freq_ten      = 1.f * num_ten  / num_repl;


  REQUIRE( std::fabs(freq_zero - expected_freq) < 30.f);
  REQUIRE( std::fabs(freq_ten - expected_freq) < 30.f);
  REQUIRE( std::fabs(freq_zero - expected_freq) < 30.f);
}

TEST_CASE( "infect_periphery" )
{
  std::cout << "testing infect periphery\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

 simulation_impl<node_2d> Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    400);


  Simulation.test_infect_periphery(1.0);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] <
                 cell_counts_before[cancer]); // after < before


  REQUIRE(cell_counts_after[infected] >  // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    76);
}

TEST_CASE( "infect_random")
{
  std::cout << "test infect random\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_2d> Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    400);

  // now we randomly infect some cells
  Simulation.test_infect_random(0.1);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] <
                 cell_counts_before[cancer]); // after < before


  REQUIRE(cell_counts_after[infected] >  // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    40);
}

TEST_CASE( "infect_center" )
{
  std::cout << "test infect center\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_2d> Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    400);

  // now we randomly infect some cells
  Simulation.test_infect_center(0.1f);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] <
                 cell_counts_before[cancer]); // after < before


  REQUIRE(cell_counts_after[infected] >  // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    40);

  size_t central_pos = Simulation.find_central_cell(cancer);
  auto ctype = Simulation.world[central_pos].get_cell_type();
  REQUIRE(ctype ==  cancer);
}

TEST_CASE( "infect_center_largest")
{
  std::cout << "test infect center\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_2d> Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  // create two squares
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  for (size_t x = 20; x < 30; ++x) {
      for (size_t y = 20; y < 30; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    500);

  Simulation.test_infect_center_largest(0.1f);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] <
                 cell_counts_before[cancer]); // after < before


  REQUIRE(cell_counts_after[infected] >  // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    40);

  auto index = Simulation.find_central_cell(infected);
  float x1 = Simulation.world[index].x_;
  float y1 = Simulation.world[index].y_;
  REQUIRE(x1 == 50);
  REQUIRE(y1 == 50);
}


TEST_CASE( "voronoi_case" )
{
  std::cout << "testing voronoi\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = true;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_2d> Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, voronoi
                                );

  REQUIRE(filler.size() ==  all_parameters.sq_num_cells *
                                   all_parameters.sq_num_cells);

  for(const auto& i : Simulation.world) {
        REQUIRE(i.neighbors.size() > 0);
  }
}

TEST_CASE( "infect_all_cancer" )
{
  // infect_all_cancer
  std::cout << "test infect all cancer\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_2d> Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    400);

  // now we randomly infect some cells
  Simulation.test_infect_all_cancer();
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] <=
                    0); // after < before


  REQUIRE(cell_counts_after[infected] >  // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    400);
}


// add_infected
TEST_CASE( "add_infected" )
{
  std::cout << "test add infected\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_2d> Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    400);

  // now we randomly infect some cells
  Simulation.add_infected(random_infection, 0.1f);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] <
                    cell_counts_before[cancer]); // after < before


  REQUIRE(cell_counts_after[infected] > // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    40);
}


// update_one_step
TEST_CASE( "update_one_step")
{
  std::cout << "test update one step\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.death_cancer = 0.0;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_2d> Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    400);

  Simulation.update_one_step();

  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();
  REQUIRE(cell_counts_after[cancer] ==
                    401);
}



// infect_long_distance
TEST_CASE( "infect_long_distance" ) {
  std::cout << "test infect long distance\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;
  all_parameters.distance_infection_upon_death = 1.0;
  all_parameters.prob_infection_upon_death = 1.f;

  simulation_impl<node_2d> Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  size_t x = 50;
  size_t y = 50;
  size_t central_pos = y * 100 + x;
  //Simulation.test_change_cell_type(central_pos, infected);
  Simulation.test_infect_long_distance(central_pos);

  static int relative_points[4][2] = { {-1, 0},
                                      {1, 0},
                                      {0, 1},
                                      {0, -1} };

  for (size_t i = 0; i < 4; ++i) {
     int x2 = x + relative_points[i][0];
     int y2 = y + relative_points[i][1];
     size_t pos = y2 * 100 + x2;
     REQUIRE(Simulation.world[pos].get_cell_type() ==
                       infected);
  }
}


// infect_long_distance
TEST_CASE( "infect_long_distance2" ) {
  std::cout << "test infect long distance\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;
  all_parameters.distance_infection_upon_death = 2.0;
  all_parameters.prob_infection_upon_death = 1.f;

  simulation_impl<node_2d> Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  int x = 50;
  int y = 50;
  int central_pos = y * 100 + x;
  //Simulation.test_change_cell_type(central_pos, infected);
  Simulation.test_infect_long_distance(central_pos);

  int relative_points[12][2] = { {-1, 0},
                                      {1, 0},
                                      {0, 1},
                                      {0, -1},
                                      {-2, 0},
                                      {2, 0},
                                      {0, -2},
                                      {0, 2},
                                      {1, 1},
                                      { 1, -1},
                                      {-1, 1},
                                      {-1, -1}};

  for (size_t i = 0; i < 12; ++i) {
     int x2 = x + relative_points[i][0];
     int y2 = y + relative_points[i][1];
     size_t pos = y2 * 100 + x2;
     REQUIRE(Simulation.world[pos].get_cell_type() ==
                       infected);
  }

  int rel_points_cancer[8][2] = { {2, 1},
                                    {2, 2},
                                    {-2, 1},
                                    {-2, 2},
                                    {2, -1},
                                    {2, -2},
                                    {-2, -1},
                                    {-2, -2}};
  for (size_t i = 0; i < 8; ++i) {
     int x2 = x + rel_points_cancer[i][0];
     int y2 = y + rel_points_cancer[i][1];
     size_t pos = y2 * 100 + x2;
     REQUIRE(Simulation.world[pos].get_cell_type() ==
                       cancer);
  }



}


TEST_CASE( "obtain_equilibrium" )
{
  std::cout << "test obtain equilibrium\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 25;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = grow;
  all_parameters.initial_number_normal_cells = 625 * 0.5;

  simulation_impl<node_2d> Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  Simulation.obtain_equilibrium(true);

  std::array<size_t, 5> ctypes = Simulation.count_cell_types();

  REQUIRE(ctypes[normal] > 0.0);
  REQUIRE(ctypes[normal] > 300);
  REQUIRE(ctypes[normal] < 10000);
  REQUIRE(ctypes[normal] < 6000);
}


TEST_CASE( "infect_second_time" )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "using random infection\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = random_infection;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 100;
  all_parameters.time_adding_virus = 100;
  all_parameters.time_adding_virus_2 = 100;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome ==  "C");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  std::string outcome2 = get_outcome(result);
  REQUIRE(outcome ==  outcome2);
}

/*
TEST_CASE( "check_grow_2" )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "Using weakened virus, expected tumor win\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.2f;
  all_parameters.death_infected = 2.0f;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 100;
  all_parameters.time_adding_virus = 200;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(result[infected] ==  0);

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  REQUIRE(result[infected] ==  0);
}
*/

TEST_CASE( "set_infection" )
{
  Param all_parameters;
  simulation_impl<node_2d> Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, regular);

  float target = 0.1f;
  Simulation.set_percent_infected(target);
  REQUIRE(Simulation.get_parameters().percent_infected ==
                    target);

  infection_routine test = random_infection;
  Simulation.set_infection_type(test);
  REQUIRE(Simulation.get_parameters().infection_type ==
                    test);
}

/*
TEST_CASE( "infect_random2" )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "using random infection\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = random_infection;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 100;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome ==  "C");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  std::string outcome2 = get_outcome(result);
  REQUIRE(outcome ==  outcome2);
}
*/
/*
TEST_CASE( "infect_periphery2" )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "using random infection\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = periphery_infection;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 100;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome ==  "C");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  std::string outcome2 = get_outcome(result);
  REQUIRE(outcome == outcome2);
}



TEST_CASE( "infect_all" )
{
  std::cout << "simulation with growth\n";
  std::cout << "testing INFECT ALL\n";
  std::cout << "using random infection\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = random_infection;
  all_parameters.percent_infected = 1.f;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 100;

  all_parameters.initial_number_normal_cells = 5000;
  all_parameters.initial_number_cancer_cells = 100;
  all_parameters.birth_cancer = 0.6f;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome ==  "A");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  std::string outcome2 = get_outcome(result);
  REQUIRE(outcome ==  outcome2);
}

TEST_CASE( "infect_center2" )
{
  std::cout << "simulation with growth\n";
  std::cout << "testing INFECT CENTER\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = center_infection;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 2;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome ==  "C");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  std::string outcome2 = get_outcome(result);
  REQUIRE(outcome ==  outcome2);
}

*/
/*
TEST_CASE( "test_converge" )
{
  std::cout << "simulation with convergence\n";
  std::cout << "Using weakened virus, expected tumor win\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.2f;
  all_parameters.death_infected = 2.0f;

  // GROW SETUP
  all_parameters.start_setup = converge;

  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 1000;
  all_parameters.time_adding_virus = 500;


  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome ==  "B");

  all_parameters.use_voronoi_grid = false;
  result = do_analysis(all_parameters);
  std::string outcome2 = get_outcome(result);
  REQUIRE(outcome ==  outcome2);
}

// TODO: add tests with long distance infection
TEST_CASE( "long_distance_infection" )
{

  std::cout << "simulation starting full\n";
  std::cout << "Using strong virus, expected virus win\n";
  std::cout << "Using long distance infection as well\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.prob_infection_upon_death = 1.f;
  all_parameters.distance_infection_upon_death = 2.0f;

  // GROW SETUP
  all_parameters.start_setup = full;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome ==  "A");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  std::string outcome2 = get_outcome(result);
  REQUIRE(outcome ==  outcome2);
}

TEST_CASE( "resistance" )
{
  std::cout << "Simulation with resistant cells\n";
  std::cout << "aiming for resistant persistence\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.prob_normal_infection = 0.01f;
  all_parameters.freq_resistant = 0.01f;

  all_parameters.start_setup = full;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
   REQUIRE(outcome ==  "D");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  std::string outcome2 = get_outcome(result);
  REQUIRE(outcome ==  outcome2);
}


TEST_CASE( "check_full_2" )
{

  std::cout << "simulation starting full\n";
  std::cout << "Using weakened virus, expected tumor win\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.2f;
  all_parameters.death_infected = 2.0f;

  // GROW SETUP
  all_parameters.start_setup = full;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome ==  "B");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  outcome = get_outcome(result);
  REQUIRE(outcome ==  "B");
}
*/
////////////////////////////////////////////////////////////////////////////////
////////////////////////  3D TESTING ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

TEST_CASE( "birth_death_3d" )
{
  std::cout << "testing birth and death in 3D\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 25;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;
 all_parameters.using_3d = true;

  simulation_impl<node_3d> Simulation(all_parameters, true);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  Simulation.t = 0.f;

  // count cell types
  std::array<size_t, 5> cells = Simulation.count_cell_types();
  for(size_t i = 0; i < 4; ++i) {
     REQUIRE(cells[i] == 0);
  }
  REQUIRE(cells[4] == all_parameters.sq_num_cells *
                              all_parameters.sq_num_cells *
                              all_parameters.sq_num_cells);


  // test update_rates

  Simulation.test_update_rates();
  for(auto c : {normal, cancer, infected, resistant}) {
   REQUIRE(Simulation.death_prob[c].get_total_sum() == 0.0);
   REQUIRE(Simulation.growth_prob[c].get_total_sum() == 0.0);
   REQUIRE(Simulation.get_rates(c) == 0.0);
  }

  Simulation.test_change_cell_type(0, normal);
  Simulation.test_update_rates();
  REQUIRE(Simulation.get_rates(1) == all_parameters.death_normal);

  Simulation.test_change_cell_type(3, cancer);
  Simulation.test_update_rates();
  REQUIRE(Simulation.get_rates(3) == all_parameters.death_cancer);

  Simulation.test_change_cell_type(0, empty);
  Simulation.test_change_cell_type(3, empty);
  Simulation.test_update_rates();

  // test change_chell_type
  // add a normal cell, then kill it

  for(auto c : {normal, cancer, infected, resistant}) {
  //  std::cout << c << "\n";
    Simulation.test_change_cell_type(5000, c);

    REQUIRE(Simulation.world[5000].get_cell_type() == c);

    REQUIRE(Simulation.death_prob[c].get_value(5000) ==
                      1.0);

    Simulation.test_change_cell_type(5000, empty);

    REQUIRE(Simulation.world[5000].get_cell_type() ==
                      empty);

    REQUIRE(Simulation.death_prob[c].get_value(5000) ==
                      0.0);
   }

  // test do_event

  std::vector< cell_type > v = {normal, empty, cancer, empty,
                                infected, empty, resistant, empty};

   for(size_t i = 0; i < 8; ++i) {
    Simulation.test_event(i);

    REQUIRE(Simulation.world[0].get_cell_type() ==
                      v[i]);

    if (v[i] != empty) {
        REQUIRE(Simulation.death_prob[v[i]].get_value(0) ==
                          1.0);
    }
   }

   // test pick event
   std::array<float, 8 > vx = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
   REQUIRE(Simulation.test_pick_event(vx, 0.0) == 0);
   for(size_t i = 0; i < 8; ++i) {
    vx[i] = 1.0f;
    REQUIRE(Simulation.test_pick_event(vx, 1.0) == i);
    vx[i] = 0.0f;
   }
}


/*
 * TODO: make 3D version
TEST_CASE( "find_central_cell_3d" )
{
  std::cout << "testing finding center cells in 3D\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;
   all_parameters.using_3d = true;

  simulation_impl<node_3d> Simulation(all_parameters, true);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  // first we try to find the central empty cell (50, 50), coordinate: 5000
  auto index = Simulation.find_central_cell(empty);
  auto x = Simulation.world[index].x_;
  auto y = Simulation.world[index].y_;
  REQUIRE(x == 49);
  REQUIRE(y == 49);

  size_t row_size = all_parameters.sq_num_cells;
  // now we add a square of cancer cells, and calculate the center
  std::vector< size_t > positions;
  for(size_t x = 10; x < 21; ++x) {
      for(size_t y = 10; y < 21; ++y) {
          size_t pos = y * row_size + x;
          Simulation.test_change_cell_type(pos, cancer);
          positions.push_back(pos);
        }
  }
  auto index2 = Simulation.find_central_cell(cancer);
  auto x2 = Simulation.world[index2].x_;
  auto y2 = Simulation.world[index2].y_;
  REQUIRE(x2 == 15);
  REQUIRE(y2 == 15);

  auto index3 = Simulation.find_central_cell(positions);
  auto x3 = Simulation.world[positions[index3]].x_;
  auto y3 = Simulation.world[positions[index3]].y_;
  REQUIRE(x2 == x3);
  REQUIRE(y2 == y3);
}
*/


TEST_CASE( "add_cells_3d" )
{
  std::cout << "testing adding cells in 3D\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;
  all_parameters.initial_number_normal_cells = 1000;
  all_parameters.initial_number_cancer_cells = 100;
 all_parameters.using_3d = true;

  simulation_impl<node_3d> Simulation(all_parameters, true);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  std::array<size_t, 5> cell_cnt = Simulation.count_cell_types();
  REQUIRE(cell_cnt[normal] == 0);

  Simulation.add_cells(normal);
  cell_cnt = Simulation.count_cell_types();
  REQUIRE(cell_cnt[normal] ==
                    all_parameters.initial_number_normal_cells);

  Simulation.add_cells(cancer);
  cell_cnt = Simulation.count_cell_types();
  REQUIRE(cell_cnt[cancer] ==
                    all_parameters.initial_number_cancer_cells);
}

TEST_CASE( "setup_types_3d" )
{
  std::cout << "testing setup types in 3D\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = full;
   all_parameters.using_3d = true;

  simulation_impl<node_3d> Simulation(all_parameters, true);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  std::array<size_t, 5> cell_cnt = Simulation.count_cell_types();

  size_t total_num_cells = all_parameters.sq_num_cells *
                           all_parameters.sq_num_cells *
                           all_parameters.sq_num_cells;

  REQUIRE(cell_cnt[normal] ==
                    total_num_cells * 0.9);
  REQUIRE(cell_cnt[cancer] ==
                    total_num_cells * 0.09);
  REQUIRE(cell_cnt[infected] ==
                    total_num_cells * 0.01);
}

TEST_CASE( "ask_infect_neighbours_3d")
{
  std::cout << "testing ask infect neighbours in 3D\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = full;
  all_parameters.distance_infection_upon_death = 1.0;
  all_parameters.prob_infection_upon_death = 100.f;
   all_parameters.using_3d = true;

  simulation_impl<node_3d> Simulation(all_parameters, true);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  size_t row_size = all_parameters.sq_num_cells;
  size_t x = 50;
  for(size_t y = 40; y < 61; ++y) {
      size_t pos = y * row_size + x;
      Simulation.test_change_cell_type(pos, cancer);
  }
  size_t y2 = 50;
  for(size_t x2 = 40; x2 < 61; ++x2) {
      size_t pos = y2 * row_size + x2;
      Simulation.test_change_cell_type(pos, cancer);
  }

  // now we have a cross of individuals
  size_t initial_pos = 50 * 100 + 50;
  Simulation.test_change_cell_type(initial_pos, infected);
  Simulation.test_ask_infect_neighbours(all_parameters.distance_infection_upon_death,
                                        initial_pos);
  for(size_t x = 40; x < 61; ++x) {
    size_t pos = 50 * row_size + x;
    int dist_x = static_cast<int>(x) - static_cast<int>(50);
    if (dist_x < 0) dist_x *= -1;
    if (dist_x == 0) {
        continue;
    }

    auto ct = Simulation.world[pos].get_cell_type();

    if (dist_x <= static_cast<int>(all_parameters.distance_infection_upon_death)) {
      REQUIRE(ct == infected);
    } else {
      REQUIRE(ct == cancer);
    }

  }
}

TEST_CASE( "random_stuff_3d" )
{
  std::cout << "testing randomizer in 3D\n";
  Param all_parameters;
   all_parameters.using_3d = true;
  all_parameters.sq_num_cells = 10;

  simulation_impl<node_3d> Simulation(all_parameters, true);

  // test random numbers
  Simulation.growth_prob[0].update_entry(10, 1.0f);
  auto x = Simulation.growth_prob[0].draw_explicit(Simulation.rndgen);
  REQUIRE(x == 10);

  Simulation.growth_prob[0].update_entry(10, 0.0f);
  Simulation.growth_prob[0].update_entry(50, 0.01f);
  auto x2 = Simulation.growth_prob[0].draw_explicit(Simulation.rndgen);
  REQUIRE(x2 == 50);

  // now we populate with equal numbers
  int total_num_cells = all_parameters.sq_num_cells * all_parameters.sq_num_cells;
  for(int i = 0; i < total_num_cells; ++i) {
    Simulation.growth_prob[0].update_entry(i, 1.0f);
  }
  // and now we draw many numbers:
  int num_zero = 0;
  int num_ten  = 0;
  int num_repl = 100000;
  for(int r = 0; r < num_repl; ++r) {
      auto xx = Simulation.growth_prob[0].draw_explicit(Simulation.rndgen);
      if(xx == 0) num_zero++;
      if(xx == 10) num_ten++;
  }

  float expected_freq = 1.f / total_num_cells;
  float freq_zero     = 1.f * num_zero / num_repl;
  float freq_ten      = 1.f * num_ten  / num_repl;


  REQUIRE(fabs(freq_zero -
                    expected_freq) < 30.f);
  REQUIRE(fabs(freq_ten -
                    expected_freq) < 30.f);
  REQUIRE(fabs(freq_zero -
                    freq_ten) < 30.f);
}

TEST_CASE( "infect_periphery_3d" )
{
  std::cout << "testing infect periphery in 3D\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;
  all_parameters.using_3d = true;

  simulation_impl<node_3d> Simulation(all_parameters, true);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          for (size_t z = 40; z < 60; ++z) {
            size_t pos = x + 100 * (y + z * 100);
            Simulation.test_change_cell_type(pos, cancer);
          }
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    8000);


  Simulation.test_infect_periphery(1.0);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] <
                 cell_counts_before[cancer]); // after < before


  REQUIRE(cell_counts_after[infected] >  // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    2168);  // 18 * 18 * 6 + 4 * 18 + 2 * 76
}

TEST_CASE( "infect_random_3d")
{
  std::cout << "test infect random in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_3d> Simulation(all_parameters, true);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          for (size_t z = 40; z < 60; ++z) {
            size_t pos = x + 100 * (y + z * 100);
            Simulation.test_change_cell_type(pos, cancer);
          }
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    8000);

  // now we randomly infect some cells
  Simulation.test_infect_random(0.1);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] <
                 cell_counts_before[cancer]); // after < before


  REQUIRE(cell_counts_after[infected] >  // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    800);
}

TEST_CASE( "infect_center_3d" )
{
  std::cout << "test infect center in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_3d> Simulation(all_parameters, true);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          for (size_t z = 40; z < 60; ++z) {
            size_t pos = x + 100 * (y + z * 100);
            Simulation.test_change_cell_type(pos, cancer);
          }
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    8000);

  // now we randomly infect some cells
  Simulation.test_infect_center(0.1f);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] <
                 cell_counts_before[cancer]); // after < before


  REQUIRE(cell_counts_after[infected] >  // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    800);

  size_t central_pos = Simulation.find_central_cell(cancer);
  auto ctype = Simulation.world[central_pos].get_cell_type();
  REQUIRE(ctype == cancer);
}

TEST_CASE( "infect_center_largest_3d")
{
  std::cout << "test infect center in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_3d> Simulation(all_parameters, true);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  // create two squares
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          for (size_t z = 40; z < 60; ++z) {
            size_t pos = x + 100 * (y + z * 100);
            Simulation.test_change_cell_type(pos, cancer);
          }
      }
  }

  for (size_t x = 20; x < 30; ++x) {
      for (size_t y = 20; y < 30; ++y) {
          for (size_t z = 20; z < 30; ++z) {
            size_t pos = x + 100 * (y + z * 100);
            Simulation.test_change_cell_type(pos, cancer);
          }
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    9000);

  Simulation.test_infect_center_largest(0.1f);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] <
                 cell_counts_before[cancer]); // after < before


  REQUIRE(cell_counts_after[infected] >  // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    800);

  auto index = Simulation.find_central_cell(infected);
  float x1 = Simulation.world[index].x_;
  float y1 = Simulation.world[index].y_;
  REQUIRE(x1 == 49.f);
  REQUIRE(y1 == 49.f);
}

TEST_CASE( "infect_all_cancer_3d" )
{
  // infect_all_cancer
  std::cout << "test infect all cancer in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_3d> Simulation(all_parameters, true);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          for (size_t z = 40; z < 60; ++z) {
            size_t pos = x + 100 * (y + z * 100);
            Simulation.test_change_cell_type(pos, cancer);
          }
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    8000);

  // now we randomly infect some cells
  Simulation.test_infect_all_cancer();
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] ==
                    0); // after < before


  REQUIRE(cell_counts_after[infected] >  // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    8000);
}


// add_infected
TEST_CASE( "add_infected_3d" )
{
  std::cout << "test add infected in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_3d> Simulation(all_parameters, true);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          for (size_t z = 40; z < 60; ++z) {
            size_t pos = x + 100 * (y + z * 100);
            Simulation.test_change_cell_type(pos, cancer);
          }
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    8000);

  // now we randomly infect some cells
  Simulation.add_infected(random_infection, 0.1f);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  REQUIRE(cell_counts_after[cancer] <
                    cell_counts_before[cancer]); // after < before


  REQUIRE(cell_counts_after[infected] >  // after > before
                 cell_counts_before[infected]);

  REQUIRE(cell_counts_after[infected] ==
                    800);
}

// update_one_step
TEST_CASE("update_one_step_3d")
{
  std::cout << "test update one step in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.death_cancer = 0.0;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation_impl<node_3d> Simulation(all_parameters, true);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          for (size_t z = 40; z < 60; ++z) {
            size_t pos = x + 100 * (y + z * 100);
            Simulation.test_change_cell_type(pos, cancer);
          }
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  REQUIRE(cell_counts_before[cancer] ==
                    8000);

  Simulation.update_one_step();

  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();
  REQUIRE(cell_counts_after[cancer] ==
                    8001);
}

/*
TEST_CASE( "obtain_equilibrium_3d ")
{
  std::cout << "test obtain equilibrium in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 25;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = grow;
  all_parameters.initial_number_normal_cells = 5000;
  all_parameters.using_3d = true;

  simulation_impl<node_3d> Simulation(all_parameters, true);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  Simulation.obtain_equilibrium();

  std::array<size_t, 5> ctypes = Simulation.count_cell_types();

  REQUIRE(ctypes[normal] > 0.0);
  REQUIRE(ctypes[normal] > 4000);
  REQUIRE(ctypes[normal] < 10000);
  REQUIRE(ctypes[normal] < 6000);
}


TEST_CASE( "infect_second_time_3d" )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "using random infection in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 25;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = random_infection;
  all_parameters.using_3d = true;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 100;
  all_parameters.time_adding_virus = 100;
  all_parameters.time_adding_virus_2 = 100;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome == "C");
}

TEST_CASE( "check_grow_2_3d" )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "Using weakened virus, expected tumor win in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 25;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.2f;
  all_parameters.death_infected = 2.0f;
  all_parameters.using_3d = true;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 100;
  all_parameters.time_adding_virus = 200;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(result[infected] == 0); // FAIL
}
*/


TEST_CASE( "set_infection_3d" )
{
  Param all_parameters;
  simulation_impl<node_3d> Simulation(all_parameters, true);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler, grid_type::regular);

  float target = 0.1f;
  Simulation.set_percent_infected(target);
  REQUIRE(Simulation.get_parameters().percent_infected ==
                    target);

  infection_routine test = random_infection;
  Simulation.set_infection_type(test);
  REQUIRE(Simulation.get_parameters().infection_type ==
                    test);
}

/*
TEST_CASE( "infect_random2_3d" )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "using random infection\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 25;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = random_infection;
  all_parameters.using_3d = true;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 100;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome == "C");
}

TEST_CASE( "infect_periphery2_3d" )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "using random infection\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 25;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = periphery_infection;
  all_parameters.using_3d = true;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 100;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome == "C");
}



TEST_CASE( "infect_all_3d" )
{
  std::cout << "simulation with growth\n";
  std::cout << "testing INFECT ALL in 3D\n";
  std::cout << "using random infection\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = random_infection;
  all_parameters.percent_infected = 1.f;
  all_parameters.using_3d = true;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 100;

  all_parameters.initial_number_normal_cells = 5000;
  all_parameters.initial_number_cancer_cells = 100;
  all_parameters.birth_cancer = 0.6f;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome == "A");
}

TEST_CASE( "infect_center2_3d" )
{
  std::cout << "simulation with growth\n";
  std::cout << "testing INFECT CENTER in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 25;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = center_infection;
  all_parameters.using_3d = true;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 2;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome == "C");
}



TEST_CASE( "test_converge_3d" )
{
  std::cout << "simulation with convergence\n";
  std::cout << "Using weakened virus, expected tumor win in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 25;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.2f;
  all_parameters.death_infected = 2.0f;
  all_parameters.using_3d = true;

  // GROW SETUP
  all_parameters.start_setup = converge;

  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 1000;
  all_parameters.time_adding_virus = 500;


  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome == "B");

}

// TODO: add tests with long distance infection
TEST_CASE( "long_distance_infection_3d" )
{

  std::cout << "simulation starting full\n";
  std::cout << "Using strong virus, expected virus win\n";
  std::cout << "Using long distance infection as well in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 25;
  all_parameters.use_voronoi_grid = false;
  all_parameters.prob_infection_upon_death = 0.5f;
  all_parameters.distance_infection_upon_death = 1.0f;
  all_parameters.using_3d = true;

  // GROW SETUP
  all_parameters.start_setup = full;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome == "A");
}

// TODO: add tests with resistant cells
TEST_CASE( "resistance_3d" )
{
  std::cout << "Simulation with resistant cells\n";
  std::cout << "aiming for resistant persistence in 3D\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 25;
  all_parameters.use_voronoi_grid = false;
  all_parameters.prob_normal_infection = 0.01f;
  all_parameters.freq_resistant = 0.01f;
  all_parameters.using_3d = true;

  all_parameters.start_setup = full;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
   REQUIRE(outcome == "D");
}


TEST_CASE( "check_full_2_3d" )
{

  std::cout << "simulation starting full\n";
  std::cout << "Using weakened virus, expected tumor win in 3D\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.2f;
  all_parameters.death_infected = 2.0f;
  all_parameters.using_3d = true;

  // GROW SETUP
  all_parameters.start_setup = full;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  REQUIRE(outcome == "B");
}
*/
