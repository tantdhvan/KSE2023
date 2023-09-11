#pragma once
#include "Constants.h"
#include <omp.h>
#include <string>
#include <vector>

using namespace std;

typedef unsigned int uint;
typedef pair<uint, uint> kpoint;
typedef vector<kpoint> kseeds;
typedef vector<double> kbins;
typedef vector<kbins> ksensors;
typedef pair<uint,double> emax;
typedef vector<emax> estimate;

class Kcommon {
public:
  Kcommon();
  ~Kcommon();

  static Kcommon *getInstance();
  unsigned randomInThread(int thread_id);

private:
  static Kcommon *instance;
  int *seed;
};
