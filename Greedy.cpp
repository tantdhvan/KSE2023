#include "Greedy.h"
#include <iostream>
#include <math.h>
#include <time.h>

Greedy::Greedy(Network *g) : Framework(g), no_queries(0) {}

Greedy::~Greedy() {}

double Greedy::get_solution(kseeds &seedsf, bool is_ds, estimate &e_max1, kpoint &fe_max, double &fe)
{
    kseeds seedst, seeds;
    no_queries = 0;
    int no_nodes = g->get_no_nodes();
    double max_inf = 0;
    vector<bool> selected(no_nodes, false);

    int e_max = -1, i_max = -1;
    double max_delta = 0;
    double max_fst = 0, max_femax = 0, max_fe2 = 0, max_fs = 0;
    int cur_i = -1;

    fe = 0;
    fe_max.first = 0;
    fe_max.second = 0;

    for (int e = 0; e < no_nodes; ++e)
    {
        max_fe2 = 0;
        cur_i = -1;
        int i_tmp = 0;
        double f_tmp = 0;
        for (int i = 0; i < Constants::K; ++i)
        {
            kseeds tmp_seeds1;
            tmp_seeds1.push_back(kpoint(e, i));
            double current_max_f = estimate_influence(tmp_seeds1);
            no_queries++;

            if (current_max_f >= f_tmp)
            {
                i_tmp = i;
                f_tmp = current_max_f;
            }

            if (current_max_f >= fe)
            {
                fe_max = kpoint(e, i);
                fe = current_max_f;
            }

            if (g->cost_matrix[e] <= Constants::BUDGET / 2)
            {
                seedst.push_back(kpoint(e, i));
                double max_fstnow = estimate_influence(seedst);
                no_queries++;
                double current_delta = max_fstnow - max_fst;
                if (current_delta >= g->cost_matrix[e] * max_fst / Constants::BUDGET)
                {
                    max_fst = max_fstnow;
                }
                else
                {
                    seedst.pop_back();
                }
            }           
        }
        e_max1[e] = s_emax(0, emax(i_tmp, f_tmp));              
    }
    double B = 0;
    vector<bool> selectedpost(no_nodes, false);

    for (int j = seedst.size() - 1; j >= 0; j--)
    {
        kpoint kp = seedst[j];
        int e = kp.first, i = kp.second;
        if (B + g->cost_matrix[e] <= Constants::BUDGET)
        {
            B += g->cost_matrix[e];
            seeds.push_back(kpoint(e, i));
        }
        else
            break;
    }

    max_fs = estimate_influence(seeds);

    no_queries++;
    seedsf.clear();
    if (max_femax <= max_fs)
    {
        seedsf = seeds;
        max_inf = max_fs;
    }
    else
    {
        seedsf.push_back(kpoint(e_max, i_max));
        max_inf = max_femax;
    }
    return max_inf;
}
double Greedy::get_solution2(kseeds &seedsf, bool is_ds)
{
    kseeds seeds;
    int no_nodes = g->get_no_nodes();
    double burned_budgets = 0.0;
    double max_inf = 0;
    int e_max = -1, i_max = -1;
    double max_femax = 0, max_fs = 0;
    vector<bool> selected(no_nodes, false);

    estimate e_max1(no_nodes);
    kpoint fe1_max;
    double fe1;
    double Gamma = Greedy::get_solution(seeds, is_ds, e_max1, fe1_max, fe1); // run Alg1

    e_max = fe1_max.first;
    i_max = fe1_max.second;
    max_femax = fe1;

    vector<double> nguong;
    int i = 0;
    while (true)
    {
        double tmp = pow((1 + Constants::EPS), i);
        if (tmp > 19 * Gamma) break;
        if(tmp>= Gamma)
            nguong.push_back(tmp);
        i++;
    }
    kseeds s;
    double max_f_s_v=0;
    for (double v : nguong)
    {
        kseeds s_v;
        double budget=0,cur_f_s_v=0;
        double theta = 2v / (5 * Constants::BUDGET);
        for(int e=0;e<no_nodes;e++)
        {
            if(budget+g->cost_matrix[e]<=Constants::BUDGET)
            {
                int cur_i=-1;
                double max_f=0;
                for(int i=0;i<Constants::K;i++)
                {
                    s_v.push_back(kpoint(e,i));
                    double f_tmp = estimate_influence(s_v);
                    no_queries++;
                    s_v.pop_back();

                    if(f_tmp>max_f)
                    {
                        cur_i=i;
                        max_f=f_tmp;
                    }
                }
                if(cur_i!=-1)
                {
                    double delta=(max_f-cur_f_s_v)/g->cost_matrix[e];
                    if(delta>=theta)
                    {
                        s_v.push_back(kpoint(e,cur_i));
                        budget+=g->cost_matrix[e];
                        cur_f_s_v=max_f;
                    }
                }
            }
        }
        if(cur_f_s_v>max_f_s_v)
        {
            seedsf=s_v;
            max_f_s_v=cur_f_s_v;
        }
    }
    cout<<"f_max: "<<max_f_s_v<<endl;
    return max(max_f_s_v,Gamma);
}
double Greedy::get_solution3(kseeds &seedsf, bool is_ds)
{
    kseeds seeds;
    no_queries = 0;
    int no_nodes = g->get_no_nodes();
    double burned_budgets = 0.0;
    double max_inf = 0;
    int e_max = -1, i_max = -1;
    double max_femax = 0, max_fs = 0;
    vector<bool> selected(no_nodes, false);

    estimate e_max1(no_nodes);
    kpoint fe1_max;
    double fe1;
    double Gamma = Greedy::get_solution(seeds, is_ds, e_max1, fe1_max, fe1);

    e_max = fe1_max.first;
    i_max = fe1_max.second;
    max_femax = fe1;

    double Theta = 10 * Gamma / (Constants::EPS * Constants::BUDGET);
    double current_f = 0.0;
    double threshold = ((1 - Constants::EPS) * Gamma) / (3 * Constants::BUDGET);
    kseeds s1;

    while (Theta >= threshold)
    {
        for (int e = 0; e < no_nodes; e++)
        {
            if (selected[e] == false && burned_budgets + g->cost_matrix[e] <= Constants::BUDGET)
            {
                double max_e_tmp = 0.0;
                int i_max_tmp = 0;
                s_emax s_e_tmp = e_max1[e];
                if (s_e_tmp.first != s1.size())
                {
                    for (int i = 0; i < Constants::K; ++i)
                    {
                        s1.push_back(kpoint(e, i));
                        double tmp = estimate_influence(s1);
                        s1.pop_back();
                        no_queries++;
                        if (tmp > max_e_tmp)
                        {
                            max_e_tmp = tmp;
                            i_max_tmp = i;
                        }
                    }
                    e_max1[e] = s_emax(s1.size(), emax(i_max_tmp, max_e_tmp)); 
                }
                else
                {
                    emax e_tmp = e_max1[e].second;
                    i_max_tmp = e_tmp.first;
                    max_e_tmp = e_tmp.second;
                }
                double delta = (max_e_tmp - current_f) / g->cost_matrix[e];
                if (delta >= Theta)
                {
                    current_f = max_e_tmp;
                    s1.push_back(kpoint(e, i_max_tmp));
                    burned_budgets += g->cost_matrix[e];
                    selected[e] = true;
                }
            }
        }
        Theta = (1 - Constants::EPS) * Theta;
    }  

    max_fs = max(max_femax, current_f);

    double l = Constants::EPS * Constants::BUDGET;
    int j = 0;
    double cost_s = 0;
    kseeds s2;
    vector<bool> check(no_nodes, false);
    double current_max_f2 = 0.0;
    while (l <= Constants::BUDGET)
    {
        int l_tmp = s2.size();
        while (j < s1.size())
        {
            kpoint kp = s1[j];
            int e = kp.first, i = kp.second;
            if ((cost_s + g->cost_matrix[e]) > l)
                break;
            cost_s += g->cost_matrix[e];
            j++;
            check[j] = true;
            s2.push_back(kp);
        }
        if (s2.size() == l_tmp)
        {
            l = (1 + Constants::EPS) * l;
            continue;
        }
        int e2 = -1, i2 = 0;
        double tmp_f2 = 0.0;
        for (int e = 0; e < no_nodes; e++)
        {
            if (check[e] == false && cost_s + g->cost_matrix[e] <= Constants::BUDGET)
            {
                double max_e_tmp = 0.0;
                int i_max_tmp = 0;
                for (int i = 0; i < Constants::K; ++i)
                {
                    s2.push_back(kpoint(e, i));
                    double tmp = estimate_influence(s2);
                    s2.pop_back();
                    no_queries++;
                    if (tmp > max_e_tmp)
                    {
                        max_e_tmp = tmp;
                        i_max_tmp = i;
                    }
                }
                if (max_e_tmp > tmp_f2)
                {
                    e2 = e;
                    i2 = i_max_tmp;
                    tmp_f2 = max_e_tmp;
                }
            }
        }
        if (e2 != -1)
        {
            if (tmp_f2 > current_max_f2)
            {
                current_max_f2 = tmp_f2;
            }
        }
        l = (1 + Constants::EPS) * l;
    }
    max_fs = max(max_fs, current_max_f2);
    return max_fs;
}
double Greedy::get_solutionT1(kseeds &seedsf)
{
    no_queries = 0;
    int no_nodes = g->get_no_nodes();
    double max_femax;
    int i_max;
    for (int e = 0; e < no_nodes; ++e)
    {
        max_femax = 0;
        i_max = 0;
        for (int i = 0; i < Constants::K; ++i)
        {
            seedsf.push_back(kpoint(e, i));
            double current_max_f = estimate_influence(seedsf);
            no_queries++;
            seedsf.pop_back();

            if (current_max_f >= max_femax)
            {
                i_max = i;
                max_femax = current_max_f;
            }
        }
        seedsf.push_back(kpoint(e, i_max));
    }
    double fmax = estimate_influence(seedsf);
    no_queries++;
    cout << fmax << endl;
    return fmax;
}
double Greedy::get_solutionT2(kseeds &seedsf)
{
    clock_t start = clock();
    kseeds Sa;
    no_queries = 0;
    double fa;
    int no_nodes = g->get_no_nodes();
    double max_f1 = 0;
    int e_max1 = 0, i_max1 = 0;
    for (int e = 0; e < no_nodes; ++e)
    {
        for (int i = 0; i < Constants::K; ++i)
        {
            if (g->cost_matrix[e] <= Constants::BUDGET)
            {
                kseeds tmp_seeds1;
                tmp_seeds1.push_back(kpoint(e, i));
                double current_max_f = estimate_influence(tmp_seeds1);
                no_queries++;
                if (current_max_f >= max_f1)
                {
                    i_max1 = i;
                    e_max1 = e;
                    max_f1 = current_max_f;
                }
            }
        }
    }
    Sa.push_back(kpoint(e_max1, i_max1));
    fa = max_f1;
    for (int e1 = 0; e1 < no_nodes - 1; ++e1)
    {
        for (int i1 = 0; i1 < Constants::K; ++i1)
        {
            for (int e2 = e1 + 1; e2 < no_nodes; ++e2)
            {
                for (int i2 = 0; i2 < Constants::K; ++i2)
                {

                    kseeds St;
                    St.push_back(kpoint(e1, i1));
                    St.push_back(kpoint(e2, i2));

                    double C_St = g->cost_matrix[e1] + g->cost_matrix[e2];
                    double f_St = estimate_influence(St);
                    no_queries++;
                    vector<bool> v(no_nodes, false);
                    v[e1] = true;
                    v[e2] = true;
                    for (int t = 0; t < no_nodes; ++t)
                    {
                        int a_max = -1, ia_max = -1;
                        double delta_max = 0, f_St_ai_max = 0;

                        for (int a = 0; a < no_nodes; ++a)
                        {
                            if (v[a] == false)
                            {
                                for (int ia = 0; ia < Constants::K; ++ia)
                                {
                                    if (((clock() - start) / CLOCKS_PER_SEC) >= 166992)
                                    {
                                        goto Timeout;
                                    }
                                    St.push_back(kpoint(a, ia));
                                    double f_St_ai = estimate_influence(St);
                                    no_queries++;
                                    St.pop_back();

                                    double delta = (f_St_ai - f_St) / g->cost_matrix[a];
                                    if (delta >= delta_max)
                                    {
                                        delta_max = delta;
                                        a_max = a;
                                        ia_max = ia;
                                        f_St_ai_max = f_St_ai;
                                    }
                                }
                            }
                        }
                        if (a_max != -1)
                        {
                            if (C_St + g->cost_matrix[a_max] <= Constants::BUDGET)
                            {
                                St.push_back(kpoint(a_max, ia_max));
                                C_St += g->cost_matrix[a_max];
                            }
                            v[a_max] = true;
                        }
                    }

                    double f_Sn = estimate_influence(St);
                    no_queries++;

                    if (f_Sn > fa)
                    {
                        Sa = St;
                        fa = f_Sn;
                    }
                }
            }
        }
    }
Timeout:
    seedsf = Sa;
    return fa;
}
int Greedy::get_no_queries()
{
    return no_queries;
}
