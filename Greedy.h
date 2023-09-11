#pragma once
#include "Framework.h"

class Greedy : public Framework
{
public:
    Greedy(Network *g);
    ~Greedy();
    double get_solution(kseeds &seedsf,bool is_ds,estimate &e_max1,kpoint &fe_max,double &fe);
    double get_solution2(kseeds &seedsf, bool is_ds);
    int get_no_queries();

private:
    int no_queries;
};
