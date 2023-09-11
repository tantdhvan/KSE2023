#include "Greedy.h"
#include <iostream>
#include <math.h>
#include <time.h>

Greedy::Greedy(Network *g) : Framework(g), no_queries(0) {}
Greedy::~Greedy() {}
double Greedy::get_solution(kseeds &seedsf,bool is_ds,estimate &e_max1,kpoint &fe_max,double &fe)
{
    kseeds seedst,seeds;
    no_queries =0;
    int no_nodes = g->get_no_nodes();
    double max_inf = 0;
    vector<bool> selected(no_nodes, false);
    int e_max = -1, i_max = -1;
    double max_delta = 0;
    double max_fst=0, max_femax=0,max_fe2=0, max_fs = 0;
    int cur_i = -1;
    fe=0;
    fe_max.first=0;
    fe_max.second=0;
    for (int e = 0; e < no_nodes; ++e)
    {
        max_fe2=0;
        cur_i=-1;
        int i_tmp=0;
        double f_tmp=0;
        for (int i = 0; i < Constants::K; ++i)
        {
            kseeds tmp_seeds1;
            tmp_seeds1.push_back(kpoint(e, i));
            double current_max_f = estimate_influence(tmp_seeds1);
            no_queries++;

            if(current_max_f>=f_tmp)
            {
                i_tmp=i;
                f_tmp=current_max_f;
            }

            if(current_max_f>=fe)
            {
                fe_max=kpoint(e,i);
                fe=current_max_f;
            }

            if (cost_matrix[e][i] > Constants::BUDGET/2)
            {
                if (current_max_f >= max_femax)
                {
                    i_max = i;
                    e_max = e;
                    max_femax = current_max_f;
                }
            }
            else
            {
                if (current_max_f >= max_fe2)
                {
                    cur_i= i;
                    max_fe2 = current_max_f;
                }
            }
        }
        e_max1[e]=emax(i_tmp,f_tmp);
        if (cur_i!=-1)
        {
            seedst.push_back(kpoint(e, cur_i));
            double max_fstnow = estimate_influence(seedst);
            no_queries++;
            double current_delta = max_fstnow - max_fst;
            if (current_delta >= cost_matrix[e][cur_i]*max_fst/Constants::BUDGET)
            {
                max_fst=max_fstnow;
            }
            else
            {
                seedst.pop_back();
            }

        }
    }
    double B=0;
    vector<bool> selectedpost(no_nodes, false);
    for(int j= seedst.size()-1; j>=0; j-- )
    {
        kpoint kp = seedst[j];
        int e = kp.first, i = kp.second;
        if (B + cost_matrix[e][i] <= Constants::BUDGET)
        {
            B +=cost_matrix[e][i];
            seeds.push_back(kpoint(e,i));
        }
        else break;
    }
    max_fs = estimate_influence(seeds);
    no_queries++;
    seedsf.clear();
    if(max_femax <= max_fs )
    {
        seedsf = seeds;
        max_inf = max_fs;
    }
    else
    {
        seedsf.push_back(kpoint(e_max,i_max));
        max_inf = max_femax;
    }
    return max_inf;
}
double Greedy::get_solution2(kseeds &seedsf,bool is_ds)
{
    kseeds seeds;
    int no_nodes = g->get_no_nodes();
    double burned_budgets = 0.0;
    double max_inf = 0;
    int e_max = -1, i_max = -1;
    double max_femax=0, max_fs = 0;
    vector<bool> selected(no_nodes, false);
    estimate e_max1(no_nodes);
    kpoint fe1_max;
    double fe1;
    double Gamma= Greedy::get_solution(seeds,is_ds,e_max1,fe1_max,fe1);
    e_max=fe1_max.first;
    i_max=fe1_max.second;
    max_femax=fe1;
    double Theta = 0;
    if(max_femax>Gamma)
    {
        seedsf.push_back(kpoint(e_max,i_max));
    }
    else
    {
        seedsf=seeds;
        max_femax=Gamma;
    }
    for(int t=0; t<=((log(10))/Constants::EPS); t++)
    {
        Theta =(5*pow((1-Constants::EPS),t)*Gamma)/Constants::BUDGET;
        if(InOutput) cout<<endl<<"Theta:"<<Theta<<endl;
        seeds.clear();
        double max_fe=0;
        for(int e=0; e<no_nodes; e++)
        {

            emax e_tmp=e_max1[e];
            int cur_i=e_tmp.first;
            double max_fi= e_tmp.second;
            seeds.push_back(kpoint(e, cur_i));
            double max_fsnow = estimate_influence(seeds);
            no_queries++;
            double current_delta = max_fsnow - max_fe;
            if ((cost_matrix[e][cur_i]+burned_budgets <= Constants::BUDGET)&&((current_delta/cost_matrix[e][cur_i]) >=Theta))
            {
                selected[e] = true;
                burned_budgets +=cost_matrix[e][cur_i];
                max_fe = max_fsnow;
            }
            else
            {
                seeds.pop_back();
            }
            if(burned_budgets>=Constants::BUDGET) break;

        }
        if( max_fe>max_femax)
        {
            seedsf=seeds;
            max_femax=max_fe;
        }
    }
    return max_femax;
}
int Greedy::get_no_queries()
{
    return no_queries;
}
