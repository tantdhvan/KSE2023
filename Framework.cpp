#include "Framework.h"
#include <math.h>
#include <iostream>

Framework::Framework(Network *g)
{
    int no_nodes = g->get_no_nodes();
    cost_matrix = Constants::cost_matrix;
    this->g = g;
    no_samples = ceil(6.75 * ((double)(g->get_no_nodes())) /
                      (Constants::EPS * Constants::EPS));
}

Framework::~Framework() {}

double Framework::estimate_test(const kseeds &seeds, uint n)
{
    double re = 0.0;
    vector<bool> a(n, false);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        a[i] = (g->sample_influence_linear_threshold(seeds) == 1);
    }
    for (bool b : a)
        if (b)
            re += 1.0;

    return re / n * g->get_no_nodes();
}

double Framework::estimate_influence(const kseeds &seeds)
{
    double re = 0.0;
    if (Constants::DATA == Social)
    {
        need_sample = ceil(((double)no_samples) / seeds.size());
        vector<bool> a(need_sample, false);
        #pragma omp parallel for
        for (int i = 0; i < need_sample; ++i)
        {
            bool tmp = (g->sample_influence_linear_threshold(seeds) == 1);

            a[i] = tmp;
        }

        for (bool b : a)
            if (b)
                re += 1.0;
        re = re / need_sample * g->get_no_nodes();

    }
    else
    {
        re = g->get_entropy(seeds);

    }
    return re;
}
