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
typedef pair<uint,double> emax;//được sử dụng để lưu cặp i,f(e) với i thuộc K sao cho f(e,i) max
typedef pair<int,emax> s_emax;//được sử dụng để lưu cặp s.size(),i,f(s+e_i)
typedef vector<s_emax> estimate;

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
