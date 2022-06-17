#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QPaintEvent>
#include <QPainter>
#include <string>
#include <sstream>
#include <chrono>
#include <iostream>

#include "../Simulation/setup.hpp"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    using_3d = false;

    sim = create_simulation(this->using_3d, all_parameters);

    basic_setup_done = false;

    ui->setupUi(this);


    ui->tab->setAutoFillBackground(true);
    ui->tab2->setAutoFillBackground(true);

    ui->line_plot->addGraph(); // normal
    ui->line_plot->addGraph(); // cancer
    ui->line_plot->addGraph(); // infected
    ui->line_plot->addGraph(); // resistant
    ui->line_plot->graph(0)->setPen(QPen(Qt::blue, 2));
    ui->line_plot->graph(1)->setPen(QPen(Qt::red, 2));
    ui->line_plot->graph(2)->setPen(QPen(Qt::green, 2));
    ui->line_plot->graph(3)->setPen(QPen(Qt::magenta, 2));

    ui->line_plot->graph(0)->setName("Normal");
    ui->line_plot->graph(1)->setName("Cancer");
    ui->line_plot->graph(2)->setName("Infected");
    ui->line_plot->graph(3)->setName("Resistant");

    QCPTextElement *fst_title = new QCPTextElement(ui->line_plot, "Number of cells");
    ui->line_plot->plotLayout()->insertRow(0);
    ui->line_plot->plotLayout()->addElement(0, 0, fst_title);
    ui->line_plot->xAxis->setLabel("Time (hours)");
    ui->line_plot->yAxis->setLabel("Number of Cells");

    ui->line_plot->legend->setVisible(true);
    QFont legendFont = font();
    legendFont.setPointSize(9);
    ui->line_plot->legend->setFont(legendFont);
    ui->line_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);

    ui->line_plot->legend->setRowSpacing(-5);
    ui->line_plot->window()->setBackgroundRole(this->window()->palette().Base);
    ui->line_plot->legend->setBrush(Qt::transparent);
    ui->line_plot->legend->setBorderPen(QPen(Qt::transparent));

    ui->t_cell_plot->addGraph();

    tcell_bars = new QCPBars(ui->t_cell_plot->xAxis,
                             ui->t_cell_plot->yAxis);

    tcell_bars->setName("Inflammation factor concentration");
    tcell_bars->setPen(QPen(Qt::blue));
    tcell_bars->setBrush(QBrush(QColor(0,0,255, static_cast<int>(0.5 * 255))));
    tcell_bars->setAntialiased(false);
    tcell_bars->setAntialiasedFill(false);
    ui->t_cell_plot->xAxis->setLabel("log10(concentration)");
    ui->t_cell_plot->yAxis->setLabel("Frequency");

    ui->progressBar->setValue(0);

    update_speed = ui->speed_slider->value();

    ui->box_infection_routine->addItem("Center");
    ui->box_infection_routine->addItem("Random");
    ui->box_infection_routine->addItem("Periphery");

    ui->box_start_setup->addItem("Full");
    ui->box_start_setup->addItem("Grow");
    ui->box_start_setup->addItem("Converge");


    ui->drpdwnbox_display->addItem("Cell types");
    ui->drpdwnbox_display->addItem("Inflammation Factor");
    ui->drpdwnbox_display->addItem("T Cell Added Death Rate");
    ui->drpdwnbox_display->addItem("Normal Growth Rate");
    ui->drpdwnbox_display->addItem("Normal Death Rate");
    ui->drpdwnbox_display->addItem("Cancer Growth Rate");
    ui->drpdwnbox_display->addItem("Cancer Death Rate");
    ui->drpdwnbox_display->addItem("Infected Growth Rate");
    ui->drpdwnbox_display->addItem("Resistant Growth Rate");
    ui->drpdwnbox_display->addItem("Dominant Growth Rate");

    ui->box_grid_type->addItem("voronoi");
    ui->box_grid_type->addItem("regular");
    ui->box_grid_type->addItem("hexagonal");

    ui->drpdwnbox_stromal_t_cell->addItem("resistant");
    ui->drpdwnbox_stromal_t_cell->addItem("sensitive");

    ui->drpdwnbox_cancer_t_cell->addItem("resistant");
    ui->drpdwnbox_cancer_t_cell->addItem("sensitive");

    ui->drpdwnbox_infected_t_cell->addItem("resistant");
    ui->drpdwnbox_infected_t_cell->addItem("sensitive");

    ui->drpdwnbox_resistant_t_cell->addItem("resistant");
    ui->drpdwnbox_resistant_t_cell->addItem("sensitive");


    is_paused = false;
    is_running = false;


    colorz.push_back(QColor(0, 0, 255));     // normal cells
    colorz.push_back(QColor(255, 0, 0));     // cancer cells
    colorz.push_back(QColor(0, 255, 0));     // infected cells
    colorz.push_back(QColor(128, 0, 128)); // purple
    colorz.push_back(QColor(0, 0, 0));       // empty cells

    ui->progressBar->setStyleSheet("QProgressBar { border-bottom-right-radius: 5px; border-bottom-left-radius: 5px; border-top-right-radius: 5px;  border-top-left-radius: 5px; border: 1px solid black; padding: 1px;background: QLinearGradient( x1: 0, y1: 0, x2: 1, y2: 0, stop: 0 #fff, stop: 1 #ddd );width: 7px;  } QProgressBar::chunk {background: QLinearGradient( x1: 0, y1: 0, x2: 1, y2: 0,stop: 0 #78d,stop: 0.4999 #46a,stop: 0.5 #45a,stop: 1 #238 );border: 1px solid black; border-bottom-right-radius: 5px; border-bottom-left-radius: 5px; border-top-right-radius: 5px;  border-top-left-radius: 5px; }");

    factor_x = 1.f;
    factor_y = 1.f;
    update_3d_gui();
    setup_simulation();
}


MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::set_resolution(int width, int height) {
    row_size = width;
    col_size = height;

    image_ = QImage(row_size, col_size, QImage::Format_RGB32);
}

void MainWindow::set_pixel(int x, int y, const QColor& col) {
    image_.setPixel(x, y, col.rgb());
}

QColor get_t_cell_color(float concentration) {
  if(concentration < 1e-5f) {
      QColor col = {0, 0, 0, 255};
      return col;
  }

  float rel_conc = concentration / 0.01; // for nicer display

  QColor col = {255, 0, 255, static_cast<int>(rel_conc * 255)};
  return col;
}

QColor get_t_cell_color_death_rate(float concentration) {
  if(concentration < 1e-5f) {
      QColor col = {0, 0, 0, 255};
      return col;
  }

  QColor col = {255, 0, 255, static_cast<int>(concentration * 255)};
  return col;
}

void MainWindow::display_regular(bool using_3d, t_cell_display display_t_cells) {
  if (!using_3d) {
    size_t line_size = static_cast<size_t>(row_size);
    size_t num_lines = static_cast<size_t>(col_size);

    for(size_t i = 0; i < num_lines; ++i) {
        QRgb* row = (QRgb*) image_.scanLine(i);

        size_t start = i * line_size;
        size_t end = start + line_size;

        for(size_t index = start; index < end; ++index) {
            size_t local_index = index - start;

            if(display_t_cells == t_cell_display::no_display) {
                row[local_index] = colorz[ sim->get_cell_type(index) ].rgb();
            }
            if(display_t_cells == t_cell_display::t_cell_rate) {
                row[local_index] = get_t_cell_color(sim->get_t_cell_concentration(index)).rgba();
              }
            if(display_t_cells == t_cell_display::added_rate) { // added death rate
               row[local_index] = get_t_cell_color_death_rate(sim->get_t_cell_rate(index) ).rgba();
            }
        }
    }
  } else {
      double frac = ui->depth_slider->value() / 100.1;
      size_t z = static_cast<size_t>(frac * row_size);

      for(size_t x = 0; x < row_size; ++x) {
          QRgb* row = (QRgb*) image_.scanLine(x);
          for(int y = 0; y < row_size; ++y) {
            int index = y + row_size * (x + z * row_size);
            //  assert(sim->world[index].z_ == z);

            if(display_t_cells == t_cell_display::no_display) {
                row[y] = colorz[ sim->get_cell_type(index) ].rgb();
            }
            if(display_t_cells == t_cell_display::t_cell_rate) {
                row[y] = get_t_cell_color(sim->get_t_cell_concentration(index)).rgba();
              }
            if(display_t_cells == t_cell_display::added_rate) { // added death rate
              // row[y] = get_t_cell_color(sim->get_added_death_rate(index) ).rgba();
                row[y] = get_t_cell_color_death_rate(sim->get_t_cell_rate(index)).rgba();
            }
          }
      }
  }
}

void MainWindow::display_voronoi(t_cell_display display_t_cells) {
  image_.fill(Qt::gray);

  QPainter painter(&image_);

  for(size_t i = 0; i < sim->num_cells; ++i) {
      QColor t_col;
      if(display_t_cells == t_cell_display::t_cell_rate) {
          t_col = get_t_cell_color(sim->get_t_cell_concentration(i));
      } else if (display_t_cells == t_cell_display::no_display) {
          t_col = colorz[sim->get_cell_type(i)];
      } else if (display_t_cells == t_cell_display::added_rate) {
          t_col = get_t_cell_color_death_rate(sim->get_t_cell_rate(i)).rgba();;
      }

      QBrush brush(t_col);
      QPainterPath path;
      path.addPolygon(polygons[i]);

      painter.fillPath(path, brush);
      painter.drawPolygon(polygons[i]);
  }
}

void MainWindow::update_image(t_cell_display display_t_cells) {

  if(grid_type == regular) {
      display_regular(using_3d, display_t_cells);
  } else {
      display_voronoi(display_t_cells);
  }

  int w = ui->q_label->width();
  int h = ui->q_label->height();

  ui->q_label->setPixmap((QPixmap::fromImage(image_)).scaled(w,h, Qt::KeepAspectRatio));
  ui->q_label->update();
}



QRgb get_color(const cell_type focal_cell_type, float rate) {

    if(rate < 1e-6f) {
        QColor col = {0, 0, 0, 255};
        return col.rgba();
    }

    if(focal_cell_type == normal) {
        QColor col = {0, 0, 255, static_cast<int>(rate * 255)};
        return col.rgba();
    }
    if(focal_cell_type == cancer) {
        QColor col = {255, 0, 0, static_cast<int>(rate * 255)};
        return col.rgba();
    }
    if(focal_cell_type == infected) {
        QColor col = {0, 255, 0, static_cast<int>(rate * 255)};
        return col.rgba();
    }
    if(focal_cell_type == resistant) {
        QColor col = {128, 0, 128, static_cast<int>(rate * 255)};
        return col.rgba();
    }

    QColor col = {0, 0, 0, 255};
    return col.rgba();
}

size_t which_max(const std::vector<float>& v) {
    size_t max_index = 0;
    float max_val = v[max_index];
    for(size_t i = 1; i < v.size(); ++i) {
        if(v[i] > max_val) {
            max_index = i;
            max_val = v[i];
        }
    }
    return max_index;
}


void MainWindow::display_regular(const binned_distribution& growth_rate,
                                 cell_type focal_cell_type,
                                 bool using_3d) {
  if (!using_3d) {
    size_t line_size = static_cast<size_t>(row_size);
    size_t num_lines = static_cast<size_t>(col_size);

    for(size_t i = 0; i < num_lines; ++i) {
        QRgb* row = (QRgb*) image_.scanLine(i);

        size_t start = i * line_size;
        size_t end = start + line_size;

        for(size_t index = start; index < end; ++index) {
            size_t local_index = index - start;

            row[local_index] = get_color(focal_cell_type, growth_rate.get_value(index));
        }
    }
  } else {
      double frac = ui->depth_slider->value() / 100.0;
      size_t z = static_cast<size_t>(frac * row_size);

      size_t num_lines = static_cast<size_t>(col_size);

      for(size_t x = 0; x < num_lines; ++x) {
          QRgb* row = (QRgb*) image_.scanLine(x);

          for(size_t y = 0; y < row_size; ++y) {
            size_t index = y + row_size * (x + z * row_size);
            //  assert(sim->world[index].z_ == z);
            row[y] = get_color(focal_cell_type, growth_rate.get_value(index));
          }
      }
  }
}

void MainWindow::display_regular(const std::array< binned_distribution, 4 > & growth_rate,
                                 bool using_3d) {

  if (!using_3d) {

    size_t line_size = static_cast<size_t>(row_size);
    size_t num_lines = static_cast<size_t>(col_size);

    for(size_t i = 0; i < num_lines; ++i) {
        QRgb* row = (QRgb*) image_.scanLine(i);

        size_t start = i * line_size;
        size_t end = start + line_size;

        for(size_t index = start; index < end; ++index) {
            size_t local_index = index - start;

            std::vector<float> probs = {0.0, 0.0, 0.0, 0.0};
            for(size_t i = 0; i < 4; ++i) {
                probs[i] = growth_rate[i].get_value(index);
            }
            cell_type focal_cell_type = static_cast<cell_type>(which_max(probs));

            row[local_index] = get_color(focal_cell_type,
                                         growth_rate[ focal_cell_type ].get_value(index));
        }
    }
  } else {
      double frac = ui->depth_slider->value() / 100.0;
      size_t z = static_cast<size_t>(frac * row_size);

      size_t num_lines = static_cast<size_t>(col_size);

      for(size_t x = 0; x < num_lines; ++x) {
          QRgb* row = (QRgb*) image_.scanLine(x);

          for(size_t y = 0; y < row_size; ++y) {
            size_t index = y + row_size * (x + z * row_size);
            //  assert(sim->world[index].z_ == z);

            std::vector<float> probs = {0.0, 0.0, 0.0, 0.0};
            for(size_t i = 0; i < 4; ++i) {
                probs[i] = growth_rate[i].get_value(index);
            }
            cell_type focal_cell_type = static_cast<cell_type>(which_max(probs));

            row[y] = get_color(focal_cell_type, growth_rate[focal_cell_type].get_value(index));
          }
      }
  }
}

void MainWindow::display_voronoi(const binned_distribution& growth_rate,
                                 cell_type focal_cell_type) {
  image_.fill(Qt::gray);

  QPainter painter(&image_);

  for(size_t i = 0; i < sim->num_cells; ++i) {

      QBrush brush(get_color(focal_cell_type,
                             growth_rate.get_value(i)));
      QPainterPath path;
      path.addPolygon(polygons[i]);

      painter.fillPath(path, brush);
      painter.drawPolygon(polygons[i]);
  }
}

void MainWindow::display_voronoi(const std::array< binned_distribution, 4 > & growth_rate) {
  image_.fill(Qt::gray);

  QPainter painter(&image_);

  for(size_t i = 0; i < sim->num_cells; ++i) {
      std::vector<float> probs = {0.0, 0.0, 0.0, 0.0};
      for(size_t k = 0; k < 4; ++k) {
          probs[k] = growth_rate[k].get_value(i);
      }
      cell_type focal_cell_type = static_cast<cell_type>(which_max(probs));

      QBrush brush(get_color(focal_cell_type,
                             probs[focal_cell_type]));
      QPainterPath path;
      path.addPolygon(polygons[i]);

      painter.fillPath(path, brush);
      painter.drawPolygon(polygons[i]);
  }
}

void MainWindow::display_regular_death_rate(const binned_distribution& death_rate,
                                            cell_type focal_cell_type,
                                            bool using_3d) {

  if (!using_3d) {
      size_t line_size = static_cast<size_t>(row_size);
      size_t num_lines = static_cast<size_t>(col_size);

      for(size_t i = 0; i < num_lines; ++i) {
          QRgb* row = (QRgb*) image_.scanLine(i);

          size_t start = i * line_size;
          size_t end = start + line_size;

          for(size_t index = start; index < end; ++index) {
              size_t local_index = index - start;

              if(focal_cell_type != cancer) {
                  row[local_index] = get_color(focal_cell_type,
                                               death_rate.get_value(index));
              } else {
                  float rate =  death_rate.get_value(index);
                  row[local_index] = get_color(focal_cell_type, rate);
              }
          }
      }
  } else {

  }
}

void MainWindow::update_image(const std::array< binned_distribution, 4 >& growth_rate) {

    cell_type focal_cell_type = normal;
    if( focal_display_type == normal_rate)    focal_cell_type = normal;
    if( focal_display_type == cancer_rate)    focal_cell_type = cancer;
    if( focal_display_type == infected_rate)  focal_cell_type = infected;
    if( focal_display_type == resistant_rate) focal_cell_type = resistant;

    if(grid_type == grid_type::regular) {
      if(focal_display_type == dominant_rate) {
          display_regular(growth_rate, using_3d);
      } else if (focal_display_type == cancer_death_rate) {
          display_regular_death_rate(growth_rate[cancer], cancer, using_3d);
      } else if (focal_display_type == normal_death_rate) {
          display_regular_death_rate(growth_rate[normal], normal, using_3d);
      } else {
          display_regular(growth_rate[focal_cell_type], focal_cell_type, using_3d);
      }
    }

    if(grid_type == grid_type::voronoi || grid_type == grid_type::hexagonal) {
        if(focal_display_type != dominant_rate) {
          display_voronoi(growth_rate[focal_cell_type], focal_cell_type);
        } else {
          display_voronoi(growth_rate);
        }
    }

    int w = ui->q_label->width();
    int h = ui->q_label->height();

    ui->q_label->setPixmap((QPixmap::fromImage(image_)).scaled(w,h, Qt::KeepAspectRatio));
    ui->q_label->update();
}

void MainWindow::update_t_cell_hist() {

    std::vector<float> vals;
    for(size_t i = 0; i < sim->num_cells; ++i) {
        if (sim->get_t_cell_concentration(i) > 0.0) {
            vals.push_back(log10(sim->get_t_cell_concentration(i)));
        }
    }
    if (vals.empty()) return;
    std::sort(vals.begin(), vals.end());
    float max_conc = vals.back();
    float min_conc = vals.front();


   int num_bins = 20;
   std::vector<size_t> hist(num_bins, 0);
   float bin_size = (max_conc - min_conc) / num_bins;
   tcell_bars->setWidth(bin_size);

   for (const auto& i : vals) {
       size_t bin = (i - min_conc) / (bin_size);
       if (bin >= hist.size()) {
        hist.back()++;
       } else {
        hist[bin]++;
       }
   }


   QVector<double> xvals(num_bins);
   QVector<double> yvals(num_bins);
 //  double max_y_val = 0.0;
   for(size_t i = 0; i < hist.size(); ++i) {
       xvals[i] = min_conc + i * bin_size;
       yvals[i] = hist[i];
 //      if(hist[i] > max_y_val) max_y_val = hist[i];
   }

   //tcell_bars->
   tcell_bars->setData(xvals, yvals);

  // ui->t_cell_plot->xAxis->setRange(-5, 0.0);

  // ui->t_cell_plot->yAxis->setRange(0, max_y_val * 1.05);

   ui->t_cell_plot->rescaleAxes();

   ui->t_cell_plot->replot();
   ui->t_cell_plot->update();
}


std::string get_string(std::string s, float v) {
    std::string output = s + " " + std::to_string(v) + "\n";
    return output;
}

void MainWindow::update_parameters(Param& p) {

   p.maximum_time = static_cast<int>(ui->box_maxtime->value());
   p.time_adding_cancer = static_cast<int>(ui->box_cancer_time->value());
   p.time_adding_virus = static_cast<int>(ui->box_virus_time->value());

   p.initial_number_cancer_cells = static_cast<size_t>(ui->box_cancer_cells->value());
   p.initial_number_normal_cells = static_cast<size_t>(ui->box_normal_cells->value());

   p.birth_normal = static_cast<float>(ui->box_birth_normal->value());
   p.death_normal = static_cast<float>(ui->box_death_normal->value());

   p.birth_cancer = static_cast<float>(ui->box_birth_cancer->value());
   p.death_cancer = static_cast<float>(ui->box_death_cancer->value());

   p.birth_infected = static_cast<float>(ui->box_birth_infected->value());
   p.death_infected = static_cast<float>(ui->box_death_infected->value());

   p.birth_cancer_resistant = static_cast<float>(ui->box_birth_cancer_resistant->value());
   p.death_cancer_resistant = static_cast<float>(ui->box_death_cancer_resistant->value());

   p.percent_infected = static_cast<float>(ui->box_percent_infected->value());
   p.prob_normal_infection = static_cast<float>(ui->box_prob_normal_infection->value());
   p.freq_resistant = static_cast<float>(ui->box_freq_resistant_cancer->value());

   p.distance_infection_upon_death = static_cast<float>(ui->box_distance_infection_death->value());
   p.prob_infection_upon_death = static_cast<float>(ui->box_prob_infection_death->value());

   p.diffusion = static_cast<float>(ui->box_diffusion->value());
   p.evaporation = static_cast<float>(ui->box_evaporation->value());
   p.t_cell_increase = static_cast<float>(ui->box_inflammation->value());
   p.t_cell_rate = static_cast<float>(ui->box_t_cell_rate->value());
   p.t_cell_density_scaler = static_cast<float>(ui->box_density_scaler->value());
   p.t_cell_inflection_point = static_cast<float>(ui->box_inflection_point->value());

   p.sq_num_cells = static_cast<size_t>(ui->box_sq_num_cells->value());

   p.infection_type = random_infection;

   auto infection_string = ui->box_infection_routine->currentText();
   if(infection_string == "Random")
       p.infection_type = random_infection;
   if(infection_string == "Center")
       p.infection_type = center_infection;
   if(infection_string == "Periphery")
       p.infection_type = periphery_infection;


    auto start_string = ui->box_start_setup->currentText();
    if(start_string == "Grow")
        p.start_setup = grow;
    if(start_string == "Full")
        p.start_setup = full;
    if(start_string == "Converge")
        p.start_setup = converge;

   auto display_string = ui->drpdwnbox_display->currentText();
   if(display_string == "Cell types")
       focal_display_type = cells;
   if(display_string == "Normal Growth Rate")
       focal_display_type = normal_rate;
   if(display_string == "Normal Death Rate")
       focal_display_type = normal_death_rate;
   if(display_string == "Cancer Growth Rate")
       focal_display_type = cancer_rate;
   if(display_string == "Cancer Death Rate")
       focal_display_type = cancer_death_rate;
   if(display_string == "Infected Growth Rate")
       focal_display_type = infected_rate;
   if(display_string == "Resistant Growth Rate")
       focal_display_type = resistant_rate;
   if(display_string == "Dominant Growth Rate")
       focal_display_type = dominant_rate;
   if(display_string == "Inflammation Factor")
       focal_display_type = inflammation_factor;
   if(display_string == "T Cell Added Death Rate")
      focal_display_type = added_death_rate;

   auto grid_string = ui->box_grid_type->currentText();
   if (grid_string == "regular") {
       grid_type = grid_type::regular;        // plotting flag
       p.use_voronoi_grid = false; // simulation flag
   }
   if (grid_string == "voronoi") {
       grid_type = grid_type::voronoi;       // plotting flag
       p.use_voronoi_grid = true; // simulation flag
   }

   if (grid_string == "hexagonal") {
       grid_type = grid_type::hexagonal; // plotting flag
       p.use_voronoi_grid = true; // simulation flag
   }

   auto t_cell_sens_stromal_string = ui->drpdwnbox_stromal_t_cell->currentText();
   if (t_cell_sens_stromal_string == "sensitive") {
       p.t_cell_sensitivity_stromal = true;
   } else {
       p.t_cell_sensitivity_stromal = false;
   }

   auto t_cell_sens_cancer_string = ui->drpdwnbox_cancer_t_cell->currentText();
   if (t_cell_sens_cancer_string == "sensitive") {
       p.t_cell_sensitivity_cancer = true;
   } else {
       p.t_cell_sensitivity_cancer = false;
   }

   auto t_cell_sens_infected_string = ui->drpdwnbox_infected_t_cell->currentText();
   if (t_cell_sens_infected_string == "sensitive") {
       p.t_cell_sensitivity_infected = true;
   } else {
       p.t_cell_sensitivity_infected = false;
   }

   auto t_cell_sens_resistant_string = ui->drpdwnbox_resistant_t_cell->currentText();
   if (t_cell_sens_resistant_string == "sensitive") {
       p.t_cell_sensitivity_resistant = true;
   } else {
       p.t_cell_sensitivity_resistant = false;
   }


   p.sq_num_pixels = static_cast<size_t>(ui->box_sq_num_pixels->value());

   if(p.use_voronoi_grid && p.sq_num_pixels <= p.sq_num_cells) {
       QMessageBox::warning(this,
                            tr("Oncolytic Virus Simulator"),
                            tr("Displaying Voronoi only works when\n"
                               "using more pixels than cells!\n"
                               "# of pixels has been increased to \n"
                               "accomodate this"));
       p.sq_num_pixels = p.sq_num_cells * 1.1f;
   }

   set_resolution(static_cast<int>(p.sq_num_pixels),
                  static_cast<int>(p.sq_num_pixels));
   return;
}

void MainWindow::update_polygons(const std::vector< std::vector< voronoi_point > >& all_edges) {

  factor_x = static_cast<float>(1.f * row_size / sim->sq_size);
  factor_y = static_cast<float>(1.f * col_size / sim->sq_size);

  polygons.clear();
  for(const auto& i : all_edges) {
      QPolygonF polygon;
      for(auto j : i) {
          polygon << QPointF(static_cast<qreal>(j.x_ * factor_x),
                             static_cast<qreal>(j.y_ * factor_y));
      }
      polygons.push_back(polygon);
  }
}

void MainWindow::setup_simulation() {
    if(is_running) {
        QMessageBox::warning(this,
                             tr("Oncolytic Virus Simulator"),
                             tr("Simulation is still running, please stop first"));
        return;
    }

    x_t.clear();
    y_n.clear();
    y_c.clear();
    y_i.clear();
    y_r.clear();

    update_parameters(all_parameters);

    sim = create_simulation(this->using_3d, all_parameters);

    //Simulation.initialize_network();
    std::vector< std::vector< voronoi_point > > all_polys;

    sim->initialize_network(all_polys,
                           grid_type);

    if(all_parameters.use_voronoi_grid == true) {
      update_polygons(all_polys);
      set_resolution(static_cast<int>(all_parameters.sq_num_pixels),
                   static_cast<int>(all_parameters.sq_num_pixels));
    }
    if(all_parameters.use_voronoi_grid == false) {
        set_resolution(static_cast<int>(all_parameters.sq_num_cells),
                     static_cast<int>(all_parameters.sq_num_cells));
    }

    sim->t = 0.0;

    ui->btn_start->setText("Start");

    update_display();

    is_paused = true;
    basic_setup_done = true;
}

void MainWindow::update_display() {

  if (focal_display_type == added_death_rate) {
      update_image(t_cell_display::added_rate);
  } else if(focal_display_type == cells) {
      update_image(t_cell_display::no_display);
  } else if (focal_display_type == inflammation_factor) {
      update_image(t_cell_display::t_cell_rate);
  } else if (focal_display_type == cancer_death_rate) {
      update_image(sim->get_death_prob());
  } else if (focal_display_type == normal_death_rate) {
      update_image(sim->get_death_prob());
  } else {
      update_image(sim->get_growth_prob());
  }

  update_plot(static_cast<double>(sim->t),
              sim->get_count_cell_types());

  if (sim->parameters.t_cell_increase > 0.0) {
      update_t_cell_hist();
  }

  QApplication::processEvents();
}



void MainWindow::obtain_equilibrium() {

  std::vector< float > densities(10, 0);
  float prev_t = sim->t;
  std::array<size_t, 5> cell_counts = sim->get_count_cell_types();
  int count = 0;
  const size_t total_num_cells = sim->num_cells;

  const static size_t range = 1 * sim->num_cells;

  int update_step = 1 + static_cast<int>((update_speed - 1) * 0.01f * range);
  int counter = 0;

  while(true) {
      sim->update_one_step();
      counter++;

      if(static_cast<int>(sim->t) - static_cast<int>(prev_t) == 10) {
           cell_counts = sim->get_count_cell_types();
           float density_normal = 1.f * cell_counts[normal] / total_num_cells;
           densities[count % 10] = cell_counts[normal];
           count++;
           if(count / 10 > 1) {
               float sum_first_half = 0.f;
               float sum_second_half = 0.f;
               for(size_t i = 0; i < 5; ++i) {
                   sum_first_half  += densities[i ];
                   sum_second_half += densities[i + 5];
               }

               std::stringstream st;
               st << sim->t << "\t" << sum_first_half * 0.2f
                            << "\t"      << sum_second_half * 0.2f << "\n";

               if(sum_first_half >= sum_second_half && density_normal > 0.4f) {
                   break;
               }
           }
           prev_t = sim->t;
      }
      if(counter % update_step == 0) {
        update_display();
      }
  }
  return;
}



void MainWindow::on_btn_start_clicked()
{
    if(is_running) return; // pre emptive return if the simulation is already running

    if(!is_paused) setup_simulation(); // only setup if it is not paused

    is_paused = false; // now everything is setup, no need to pause.
    ui->btn_start->setText("Start");

    is_running = true;
    int counter = 0;

    if (all_parameters.start_setup == converge) {
        obtain_equilibrium();
        sim->t = 0.f;
        sim->set_start_setup(grow);
        x_t.clear();
        y_n.clear();
        y_c.clear();
        y_i.clear();
        y_r.clear();
        sim->update_grow_params();
    }

    while(sim->t < all_parameters.maximum_time) {
        sim->update_one_step();
        counter++;

        int progress = static_cast<int>(100.f * sim->t / all_parameters.maximum_time);
        ui->progressBar->setValue(progress);

        // update speed is in 1 - 100
        const static size_t range = 1 * sim->num_cells;

        // speed-1 * 0.01f to normalise value in [0, 1] (original value is [1, 100]
        // then times range, maximum speed is thus total number of cells.
        int update_step = 1 + static_cast<int>((update_speed - 1) * 0.01f * range);
        if (update_speed > 90) {
            // create extreme speed.
          update_step = 1 + static_cast<int>((update_speed - 1) * 0.1f * range);
        }
        if (update_speed < 5) {
            update_step = update_speed;
        }

        if(counter % update_step == 0) {
            update_display();
        }
        if(!is_running) break;
    }
    is_running = false;
}

void MainWindow::update_plot(double t,
                             const std::array<size_t, 5>& cell_numbers) {
    x_t.append(t);
    y_n.append(cell_numbers[0]);
    y_c.append(cell_numbers[1]);
    y_i.append(cell_numbers[2]);
    y_r.append(cell_numbers[3]);

    ui->line_plot->graph(0)->data()->clear();
    ui->line_plot->graph(0)->setData(x_t, y_n);

    ui->line_plot->graph(1)->data()->clear();
    ui->line_plot->graph(1)->setData(x_t, y_c);

    ui->line_plot->graph(2)->data()->clear();
    ui->line_plot->graph(2)->setData(x_t, y_i);

    ui->line_plot->graph(3)->data()->clear();
    ui->line_plot->graph(3)->setData(x_t, y_r);


    ui->line_plot->rescaleAxes();
    ui->line_plot->replot();
}

void MainWindow::on_btn_stop_clicked() {
    is_running = false;
    is_paused = true;
    ui->btn_start->setText("Resume");
}

void MainWindow::on_speed_slider_actionTriggered(int) {
   update_speed = ui->speed_slider->value();
}

void MainWindow::on_drpdwnbox_display_activated(int)
{
    auto display_string = ui->drpdwnbox_display->currentText();
    if(display_string == "Cell types")
        focal_display_type = cells;
    if(display_string == "Normal Growth Rate")
        focal_display_type = normal_rate;
    if(display_string == "Normal Death Rate")
        focal_display_type = normal_death_rate;
    if(display_string == "Cancer Growth Rate")
        focal_display_type = cancer_rate;
    if(display_string == "Cancer Death Rate")
        focal_display_type = cancer_death_rate;
    if(display_string == "Infected Growth Rate")
        focal_display_type = infected_rate;
    if(display_string == "Resistant Growth Rate")
        focal_display_type = resistant_rate;
    if(display_string == "Dominant Growth Rate")
        focal_display_type = dominant_rate;
    if(display_string == "Inflammation Factor")
        focal_display_type = inflammation_factor;
    if(display_string == "T Cell Added Death Rate")
       focal_display_type = added_death_rate;

    if (!is_running) {
       update_display();
       QApplication::processEvents();
    }
}

void MainWindow::on_btn_setup_clicked()
{
  ui->progressBar->setValue(0);
  setup_simulation();
}

void MainWindow::on_btn_add_virus_clicked()
{
   auto infection_string = ui->box_infection_routine->currentText();
   if(infection_string == "Random")
       sim->set_infection_type(random_infection);
   if(infection_string == "Center")
       sim->set_infection_type(center_infection);
   if(infection_string == "Periphery")
       sim->set_infection_type(periphery_infection);

   sim->set_percent_infected(static_cast<float>(ui->box_percent_infected->value()));

   sim->add_infected(sim->get_infection_type(),
                     sim->get_percent_infected());
}

void MainWindow::update_3d_gui() {
   ui->depth_slider->setVisible(using_3d); // depth slider
   ui->label_30->setVisible(using_3d);     // Focal Layer
   ui->label_31->setVisible(using_3d);     // Front
   ui->label_32->setVisible(using_3d);     // Middle
   ui->label_33->setVisible(using_3d);     // Back
}

void MainWindow::on_btn_use_3d_clicked()
{
  using_3d = !using_3d;
  update_3d_gui();
  if (using_3d) {
      ui->btn_use_3d->setText("Switch to 2D");
      // ui->label_using_3d->setText("USING 3D");

      if (grid_type != grid_type::regular) {
          QMessageBox::warning(this,
                               tr("Oncolytic Virus Simulator"),
                               tr("Using 3-Dimensions only works\n"
                                  "with a regular grid, now setting\n"
                                  "grid to regular!"));
      }
      grid_type = grid_type::regular;
      all_parameters.use_voronoi_grid = false;
      ui->box_grid_type->setCurrentText("regular");

  } else {
      ui->btn_use_3d->setText("Switch to 3D");
   //   ui->label_using_3d->setText("using 2D, NOT using 3D");
  }


  ui->progressBar->setValue(0);
  is_running = false;
  is_paused = true;
  setup_simulation();
}

// boxes!
void MainWindow::on_box_death_cancer_resistant_valueChanged()
{
    sim->parameters.death_cancer_resistant = ui->box_death_cancer_resistant->value();
}

void MainWindow::on_box_birth_infected_valueChanged()
{
  sim->parameters.birth_infected = ui->box_birth_infected->value();
}

void MainWindow::on_box_death_infected_valueChanged()
{
    sim->parameters.death_infected = ui->box_death_infected->value();
}

void MainWindow::on_box_birth_cancer_resistant_valueChanged()
{
  sim->parameters.birth_cancer_resistant = ui->box_birth_cancer_resistant->value();
}

void MainWindow::on_box_grid_type_currentIndexChanged()
{
  if (basic_setup_done) { // we don't want this to trigger prematurely
    auto grid_string = ui->box_grid_type->currentText();
    if (grid_string == "regular") {
        grid_type = grid_type::regular;        // plotting flag
        sim->parameters.use_voronoi_grid = false; // simulation flag
    }
    if (grid_string == "voronoi") {
        grid_type = grid_type::voronoi;       // plotting flag
        sim->parameters.use_voronoi_grid = true; // simulation flag
    }

    if (grid_string == "hexagonal") {
        grid_type = grid_type::hexagonal; // plotting flag
        sim->parameters.use_voronoi_grid = true; // simulation flag
    }
    is_running = false;
    is_paused = true;
    setup_simulation();
  }
}

void MainWindow::on_box_sq_num_cells_valueChanged()
{
    if (ui->box_sq_num_cells->value() > 10) {
      all_parameters.sq_num_cells = ui->box_sq_num_cells->value();
      is_running = false;
      is_paused = true;
      setup_simulation();
    }
}

void MainWindow::on_box_birth_normal_valueChanged()
{
    sim->parameters.birth_normal = ui->box_birth_normal->value();
}

void MainWindow::on_box_death_normal_valueChanged()
{
    sim->parameters.death_normal = ui->box_death_normal->value();
}

void MainWindow::on_box_birth_cancer_valueChanged()
{
    sim->parameters.birth_cancer = ui->box_death_cancer->value();
}

void MainWindow::on_box_death_cancer_valueChanged()
{
    sim->parameters.death_cancer = ui->box_death_cancer->value();
}

void MainWindow::on_box_maxtime_valueChanged()
{
    sim->parameters.maximum_time = ui->box_maxtime->value();
}

void MainWindow::on_box_cancer_time_valueChanged()
{
    sim->parameters.time_adding_cancer = ui->box_cancer_time->value();
}

void MainWindow::on_box_virus_time_valueChanged()
{
  sim->parameters.time_adding_virus = ui->box_virus_time->value();
}

void MainWindow::on_box_prob_normal_infection_valueChanged()
{
  sim->parameters.prob_normal_infection = ui->box_prob_normal_infection->value();
}

void MainWindow::on_box_freq_resistant_cancer_valueChanged()
{
    sim->parameters.freq_resistant = ui->box_freq_resistant_cancer->value();
}

void MainWindow::on_box_percent_infected_valueChanged()
{
    sim->parameters.percent_infected = ui->box_percent_infected->value();
}

void MainWindow::on_box_sq_num_pixels_valueChanged()
{
    sim->parameters.sq_num_pixels = ui->box_sq_num_pixels->value();
}

void MainWindow::on_box_normal_cells_valueChanged()
{
    sim->parameters.initial_number_normal_cells = ui->box_normal_cells->value();
}

void MainWindow::on_box_cancer_cells_valueChanged()
{
    sim->parameters.initial_number_cancer_cells = ui->box_cancer_cells->value();
}

void MainWindow::on_box_distance_infection_death_valueChanged()
{
    sim->parameters.distance_infection_upon_death = ui->box_distance_infection_death->value();
}

void MainWindow::on_box_prob_infection_death_valueChanged()
{
    sim->parameters.prob_infection_upon_death = ui->box_prob_infection_death->value();
}

void MainWindow::on_box_infection_routine_currentIndexChanged()
{
    auto inf_type = ui->box_infection_routine->currentText();
    if (inf_type == "Center") {
         sim->parameters.infection_type = infection_routine::center_infection;
    }
    if (inf_type == "Random") {
        sim->parameters.infection_type = infection_routine::random_infection;
    }
    if (inf_type == "Periphery") {
        sim->parameters.infection_type = infection_routine::periphery_infection;
    }

    if (sim->t > sim->parameters.time_adding_virus) {
        is_running = false;
        is_paused = true;
        setup_simulation();
    }
}

void MainWindow::on_depth_slider_actionTriggered()
{
    update_display();
}

void MainWindow::on_depth_slider_sliderMoved()
{
    update_display();
}

void MainWindow::on_box_diffusion_valueChanged()
{
   sim->parameters.diffusion = ui->box_diffusion->value();
}

void MainWindow::on_box_evaporation_valueChanged()
{
    sim->parameters.evaporation = ui->box_evaporation->value();
}

void MainWindow::on_box_inflammation_valueChanged()
{
    sim->parameters.t_cell_increase = ui->box_inflammation->value();
}

void MainWindow::on_box_t_cell_rate_valueChanged()
{
    sim->parameters.t_cell_rate = ui->box_t_cell_rate->value();
}

void MainWindow::on_box_density_scaler_valueChanged()
{
    sim->parameters.t_cell_density_scaler = ui->box_density_scaler->value();
}

void MainWindow::on_box_inflection_point_valueChanged()
{
    sim->parameters.t_cell_inflection_point = ui->box_inflection_point->value();
}

void MainWindow::on_box_start_setup_activated()
{
  auto start_string = ui->box_start_setup->currentText();
      if(start_string == "Grow")
          sim->parameters.start_setup = grow;
      if(start_string == "Full")
          sim->parameters.start_setup = full;
      if(start_string == "Converge")
          sim->parameters.start_setup = converge;

      is_running = false;
      is_paused = true;
      setup_simulation();
}

void MainWindow::on_drpdwnbox_stromal_t_cell_activated(int index)
{
    auto t_cell_sens_stromal_string = ui->drpdwnbox_stromal_t_cell->currentText();
    if (t_cell_sens_stromal_string == "sensitive") {
        sim->parameters.t_cell_sensitivity_stromal = true;
    } else {
        sim->parameters.t_cell_sensitivity_stromal = false;
    }
}


void MainWindow::on_drpdwnbox_cancer_t_cell_activated(int index)
{
    auto t_cell_sens_cancer_string = ui->drpdwnbox_cancer_t_cell->currentText();
    if (t_cell_sens_cancer_string == "sensitive") {
        sim->parameters.t_cell_sensitivity_cancer = true;
    } else {
        sim->parameters.t_cell_sensitivity_cancer = false;
    }
}


void MainWindow::on_drpdwnbox_infected_t_cell_activated(int index)
{
    auto t_cell_sens_infected_string = ui->drpdwnbox_infected_t_cell->currentText();
    if (t_cell_sens_infected_string == "sensitive") {
        sim->parameters.t_cell_sensitivity_infected = true;
    } else {
        sim->parameters.t_cell_sensitivity_infected = false;
    }
}


void MainWindow::on_drpdwnbox_resistant_t_cell_activated(int index)
{
    auto t_cell_sens_resistant_string = ui->drpdwnbox_resistant_t_cell->currentText();
    if (t_cell_sens_resistant_string == "sensitive") {
        sim->parameters.t_cell_sensitivity_resistant = true;
    } else {
        sim->parameters.t_cell_sensitivity_resistant = false;
    }
}

