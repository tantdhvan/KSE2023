#pragma once
#include "Network.h"
#include <vector>

class Framework {
public:
  Framework(Network *g);
  ~Framework();
  virtual double get_solution(kseeds &seedsf,bool is_ds,estimate &e_max1,kpoint &fe_max,double &fe) = 0;
  double estimate_influence(const kseeds &seeds);
  double estimate_test(const kseeds &seeds, uint n);

protected:
  Network *g;
  uint no_samples;
  vector<vector<double>> cost_matrix;
};
