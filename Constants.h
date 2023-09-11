#pragma once
#include <vector>

enum DATATYPE { Social, Sensor };
enum ALG { Adaptive, eAppro, ream, SGr, SampleRstream, Tang, Tang2 };

class Constants {
public:
  Constants();
  ~Constants();

  static ALG ALGORITHM;
  static DATATYPE DATA;
  static int NUM_THREAD;
  static double BUDGET;
  static int K;
  static int M;
  static double EPS;
  static double EPS_TAG;
  static double ALPHA;
  static double BETA;
  static double RHO;
  static double DELTA;         // for streaming algorithm
  static int NO_DENOISE_STEPS; // steps for eps' = 0 reach epsilon
  static std::vector<std::vector<double>> cost_matrix;
};
