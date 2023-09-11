#include "Streaming.h"
#include <algorithm>
#include <iostream>
#include <math.h>

Streaming::Streaming(Network *g) : Framework(g)
{
  int no_nodes = g->get_no_nodes();
  for (int i = 0; i < no_nodes; ++i)
  {
    node_sequence.push_back(i);
  }
  random_shuffle(node_sequence.begin(), node_sequence.end());
  no_queries = 0;
  max_solution = 0.0;
}

Streaming::~Streaming() {}

int Streaming::select_element_ds(int j, uint e, int step)
{
  uint i_max = 0;
  double max_inf = 0.0;

  for (int i = 0; i < Constants::K; ++i)
  {
    kseeds seeds_tmp = sub_seeds[j][step];
    seeds_tmp.push_back(kpoint(e, i));
    double current_inf = estimate_influence(seeds_tmp);
    ++no_queries;
    if (max_inf < current_inf)
    {
      max_inf = current_inf;
      i_max = i;
    }
  }

  /*
   * No need to consider denoise step in deterministic
   * double eps = Constants::EPS / (Constants::NO_DENOISE_STEPS > 1 ?
   * Constants::NO_DENOISE_STEPS - 1 : Constants::NO_DENOISE_STEPS) * step;
   */
  if (i_max >= 0 &&
      (max_inf / (sub_seeds_cost[j][step] + g->cost_matrix[e])) >=
          (Constants::ALPHA * thresholds[j] / Constants::BUDGET))
  {
    max_solution = max(max_solution, max_inf);
    //cout<<"max_solution: "<<max_solution<<endl;
    return i_max;
  }
  else
  {
    return -1;
  }
}
int Streaming::select_element_rs(int j, uint e, int step)
{
  vector<double> temp_inf(Constants::K), p(Constants::K);
  uint J_size = 0;

  for (int i = 0; i < Constants::K; ++i)
  {
    kseeds tmp_seeds = sub_seeds[j][step];
    tmp_seeds.push_back(kpoint(e, i));
    temp_inf[i] = estimate_influence(tmp_seeds);
    ++no_queries;

    bool critical_cond =
        sub_seeds_cost[j][step] + g->cost_matrix[e] <= Constants::BUDGET &&
        ((temp_inf[i] - influences[j][step]) / g->cost_matrix[e] >=
         (Constants::ALPHA * thresholds[j]) / Constants::BUDGET);

    if (critical_cond)
    {
      p[i] = (temp_inf[i] - influences[j][step]) / g->cost_matrix[e];
      if (p[i] > 0)
        ++J_size;
    }
  }

  if (J_size == 0)
  {
    return -1;
  }
  else if (J_size == 1)
  {
    for (int i = 0; i < Constants::K; ++i)
    {
      if (p[i] > 0)
      {
        influences[j][step] = temp_inf[i];
        if (influences[j][step] > max_solution)
          max_solution = influences[j][step];
          //cout<<"max_solution:"<<max_solution<<endl;
        return i;
      }
    }
  }
  else
  {
    double T = 0.0;
    double sum = 0.0;

    for (int i = 0; i < Constants::K; ++i)
    {
      p[i] = pow(p[i], J_size - 1);
      T += p[i];
    }

    double random = (double)(rand() % 10000) / 10000 * T;
    for (int i = 0; i < Constants::K; ++i)
    {
      if (sum <= random && random < sum + p[i])
      {
        influences[j][step] = temp_inf[i];
        max_solution = max(max_solution, influences[j][step]);
        //cout<<"max_solution:"<<max_solution<<endl;
        return i;
      }
      sum += p[i];
    }
  }

  return -1;
}
double Streaming::get_solution(bool is_dstream)
{
  int j_min, j_max = 0, p_max = 0;
  int c_passed = 0;
  double max_gain = 0.0;
  double log_delta = log(1 + Constants::EPS_TAG);

  for (uint e : node_sequence)
  {
    ++c_passed;
    for (int i = 0; i < Constants::K; ++i)
    {
      kseeds single_seed = {kpoint(e, i)};
      double current_inf = estimate_influence(single_seed);
      ++no_queries;

      if (max_gain < current_inf)
      {
        max_gain = current_inf;
        j_min = ceil(log(max_gain / log_delta));

        if (is_dstream)
        {
          /* upper bound if ds is running */
          j_max = floor(log(max_gain * Constants::BUDGET) / log_delta);
        }
        else
        {
          /* upper bound if rs is running */
          j_max = floor(log(max_gain * Constants::BUDGET) / log_delta);
        }

        j_min = max(j_min, 0); /* could < 0 */
        j_max = max(j_max, 0); /* could < 0 */

        /* Resize subseeds if having new j variable */
        if (j_max + 1 > sub_seeds.size())
        {
          sub_seeds.resize(j_max + 1,
                           vector<kseeds>(Constants::NO_DENOISE_STEPS));
          sub_seeds_cost.resize(
              j_max + 1, vector<double>(Constants::NO_DENOISE_STEPS, 0.0));
          influences.resize(j_max + 1,
                            vector<double>(Constants::NO_DENOISE_STEPS, 0.0));
          for (int j = thresholds.size(); j < j_max + 1; ++j)
          {
            thresholds.push_back(pow(1 + Constants::EPS_TAG, j));
          }
        }
      }
    }

    /* Special case for Algorithm 5 */
    vector<int> rhos(j_max + 1);
    if (Constants::NO_DENOISE_STEPS > 1)
    {
      for (int j = j_min; j <= j_max; ++j)
      {
        int p = floor(
            (4 * log(1 / Constants::DELTA)) /
            (pow(Constants::EPS_TAG, 2) * Constants::ALPHA * thresholds[j]));
        cout << "/* Special case for Algorithm 5 */"
             << " " << p << endl;
        rhos[j] = max(0, p);
        if (p + 1 > sub_seeds_cost[j].size())
        {
          sub_seeds[j].resize(p + 1, vector<kpoint>(0));
          sub_seeds_cost[j].resize(p + 1, 0.0);
          influences[j].resize(p + 1, 0.0);
        }
      }
    }

    if (sub_seeds.size() > 0)
    {
      for (int j = j_min; j <= j_max; ++j)
      {
        int denoise_step = max(1, rhos[j]);
        for (int step = 0; step < denoise_step; ++step)
        {
          if (sub_seeds_cost[j][step] + g->cost_matrix[e] < Constants::BUDGET)
          {
            int i;
            if (is_dstream)
              i = select_element_ds(j, e, step);
            else
              i = select_element_rs(j, e, step);
            //cout<<"max_solution_new:"<<max_solution<<endl;
            if (i != -1)
            {
              sub_seeds[j][step].push_back(kpoint(e, i));
              sub_seeds_cost[j][step] += g->cost_matrix[e];
            }
          }
        }
      }
    }
  }
  cout<<"max_solution_new:"<<max_solution<<endl;
  return max_solution;
}

uint Streaming::get_no_queries() { return no_queries; }
