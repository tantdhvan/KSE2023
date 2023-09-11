#pragma once
#include "Framework.h"

class Streaming : public Framework {
public:
  Streaming(Network *g);
  ~Streaming();
  double get_solution(bool is_ds = true);
  int select_element_ds(int j, uint e, int step);
  int select_element_rs(int j, uint e, int step);
  uint get_no_queries();

protected:
  double max_solution;
  vector<vector<kseeds>> sub_seeds;
  vector<vector<double>> sub_seeds_cost;
  vector<double> thresholds;
  vector<vector<double>> influences;
  vector<uint> node_sequence;
  uint no_queries;
};
