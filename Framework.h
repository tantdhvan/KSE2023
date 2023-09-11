#pragma once
#include "Network.h"
#include <vector>

class Framework {
public:
  Framework(Network *g);
  ~Framework();
  double estimate_influence(const kseeds &seeds);
  double estimate_test(const kseeds &seeds, uint n);
protected:
  Network *g;
  uint no_samples;
};
