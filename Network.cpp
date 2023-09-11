#include "Network.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <numeric>
#include <queue>
#include <sstream>
#include <random>

using namespace std;

Network::Network()
{
    common_instance = Kcommon::getInstance();
}

Network::~Network() {}

int Network::get_no_nodes()
{
    return number_of_nodes;
}

int Network::get_out_degree(uint n)
{
    return out_neighbors[n].size();
}
void Network::init_alpha()
{
    uniform_real_distribution<double> unidist(0, 1);
    alpha.assign(number_of_nodes, 0.0);
    mt19937 gen(0); // same sequence each time
    for (int u = 0; u < number_of_nodes; ++u)
    {
        alpha[u] = unidist(gen);
    }
}
bool Network::read_network_from_file(int no_nodes, string file, bool is_directed)
{
    clear();

    number_of_nodes = no_nodes;
    out_neighbors = vector<vector<uint>>(no_nodes);
    in_neighbors = vector<vector<uint>>(no_nodes);

    // initiate preferences
    preferences = vector<vector<uint>>(no_nodes);
    vector<uint> ini(Constants::K);

    for (int i = 0; i < Constants::K; ++i)
    {
        ini[i] = i;
        probabilities.push_back(1.0 / (2 * Constants::K) * (1 + i)); // probability for each level of reference
        // is 1/2k, 2/2k, 3/2k, ... k/2k
    }

    for (int i = 0; i < no_nodes; ++i)
    {
        preferences[i] = ini;
        // random shuffle
        random_shuffle(preferences[i].begin(), preferences[i].end());
    }

    this->is_directed = is_directed;
    ifstream is(file);
    is.seekg(0, is.end);
    long bufSize = is.tellg();
    is.seekg(0, is.beg);
    int item = 0;

    char *buffer = new char[bufSize];

    is.read(buffer, bufSize);
    is.close();

    std::string::size_type sz = 0;
    long sp = 0;
    uint start_id, end_id;
    bool is_start = true;
    int zero_id = -1;
    uint id = 0;
    uint s_id, e_id; // used to stored ordered id of startId and endId
    uint edge_id = 0;

    while (sp < bufSize)
    {
        char c = buffer[sp];
        item = item * 10 + c - 48;
        sp++;
        if (sp == bufSize ||
            (sp < bufSize && (buffer[sp] < 48 || buffer[sp] > 57)))
        {
            while (sp < bufSize && (buffer[sp] < 48 || buffer[sp] > 57))
            {
                sp++;
            }

            if (is_start)
            {
                start_id = item;
                if (start_id >= no_nodes)
                    return false;
                is_start = false;
            }
            else
            {
                end_id = item;
                if (end_id >= no_nodes)
                    return false;
                is_start = true;

                s_id = map_node_id[start_id];
                if (!s_id && start_id != zero_id) // if start node has not exist
                {
                    map_node_id[start_id] = id;
                    s_id = id;
                    if (id == 0)
                        zero_id = start_id;
                    id++;
                }

                e_id = map_node_id[end_id];
                if (!e_id && end_id != zero_id) // if end node has not exist
                {
                    map_node_id[end_id] = id;
                    e_id = id;
                    id++;
                }

                // add out neighbor
                out_neighbors[s_id].push_back(e_id);
                in_neighbors[e_id].push_back(s_id);

                if (!is_directed)
                {
                    out_neighbors[e_id].push_back(s_id);
                    in_neighbors[s_id].push_back(e_id);
                }
            }

            item = 0;
        }
    }

    delete[] buffer;
    return true;
}

bool Network::read_network_from_file2(string file, bool is_directed)
{
    clear();
    int no_edges, no_nodes;
    ifstream infile(file);
    infile >> no_nodes >> no_edges;

    number_of_nodes = no_nodes;
    my_neighbors = vector<vector<pair<int, double>>>(no_nodes);
    my_weight_k = vector<vector<double>>(no_nodes);
    int s_node, e_node;
    double w;
    for (int i = 0; i < no_edges; i++)
    {
        infile >> s_node >> e_node >> w;
        my_neighbors[s_node].push_back(make_pair(e_node, w));
        my_neighbors[e_node].push_back(make_pair(s_node, w));
    }
    double k1, k2, k3;
    for (int i = 0; i < no_nodes; i++)
    {
        infile >> k1 >> k2 >> k3;
        my_weight_k[i].push_back(k1);
        my_weight_k[i].push_back(k2);
        my_weight_k[i].push_back(k3);
    }
    double c;
    cost_matrix=vector<double>(no_nodes);
    for (int i = 0; i < no_nodes; i++)
    {
        infile >> c;
        cost_matrix[i]=c;
    }
    alpha=vector<double>(no_nodes);
    double al;
    for (int i = 0; i < no_nodes; i++)
    {
        infile >> al;
        alpha[i]=al;
    }
    return true;
}

void Network::generate_random_network(int no_nodes, double p, bool is_directed)
{
    clear();
    number_of_nodes = no_nodes;
    this->is_directed = is_directed;
    out_neighbors = vector<vector<uint>>(no_nodes);
    in_neighbors = vector<vector<uint>>(no_nodes);

    // initiate preferences
    preferences = vector<vector<uint>>(no_nodes);
    vector<uint> ini(Constants::K);

    for (int i = 0; i < Constants::K; ++i)
    {
        ini[i] = i;
        probabilities.push_back(1.0 / (2 * Constants::K) *
                                (1 + i)); // probability for each level of reference
        // is 1/2k, 2/2k, 3/2k, ... k/2k
    }

    for (int i = 0; i < no_nodes; ++i)
    {
        preferences[i] = ini;
        // random shuffle
        random_shuffle(preferences[i].begin(), preferences[i].end());
    }

#pragma omp parallel for
    for (int i = 0; i < number_of_nodes; ++i)
    {
        for (int j = is_directed ? 0 : i + 1; j < number_of_nodes; ++j)
        {
            if ((is_directed && i != j) || !is_directed)
            {
                int r = common_instance->randomInThread(omp_get_thread_num()) % 1000;
                if (r < p * 1000)
                {
#pragma omp critical
                    {
                        out_neighbors[i].push_back(j);
                        in_neighbors[j].push_back(i);

                        if (!is_directed)
                        {
                            out_neighbors[j].push_back(i);
                            in_neighbors[i].push_back(j);
                        }
                    }
                }
            }
        }
    }
}

uint Network::sample_influence(const kseeds &seeds)
{
    uint re = 0;
    queue<pair<uint, int>> q;
    vector<bool> visited(number_of_nodes, false);
    for (pair<uint, int> seed : seeds)
    {
        q.push(seed);
        visited[seed.first] = true;
    }

    while (!q.empty())
    {
        pair<uint, int> p = q.front();
        ++re;
        q.pop();
        double prob = 0;    // probability to influence an out-neighbor
        if (p.second != -1) // this pair is from seeds
        {
            prob = ((double)preferences[p.first][p.second]) /
                   out_neighbors[p.first].size();
        }
        else
        {
            prob = 1.0 / out_neighbors[p.first].size();
        }

        for (uint nei : out_neighbors[p.first])
        {
            if (!visited[nei])
            {
                double r =
                    ((double)(common_instance->randomInThread(omp_get_thread_num()) %
                              10000)) /
                    10000.0;
                if (r <= prob)
                {
                    q.push(pair<uint, int>(nei, -1));
                    visited[nei] = true;
                }
            }
        }
    }

    return re;
}

uint Network::sample_influence_reverse(const kseeds &seeds)
{
    vector<int> pref(number_of_nodes, -1);
    for (pair<uint, int> p : seeds)
    {
        pref[p.first] = p.second;
    }

    int r_node =
        common_instance->randomInThread(omp_get_thread_num()) % number_of_nodes;
    if (pref[r_node] != -1) // if random node is of seeds
    {
        return 1;
    }
    else
    {
        queue<uint> q;
        vector<bool> visited(number_of_nodes, false);
        q.push(r_node);
        visited[r_node] = true;
        while (!q.empty())
        {
            uint u = q.front();
            q.pop();
            for (uint v : in_neighbors[u])
            {
                if (!visited[v])
                {
                    double prob = 0.0;
                    if (pref[v] != -1) // this neighbor is of seeds
                    {
                        prob = ((double)preferences[v][pref[v]] + 1.0) /
                               out_neighbors[v].size();
                        // prob = probabilities[preferences[v][pref[v]]];
                    }
                    else
                    {
                        prob = 1.0 / out_neighbors[v].size();
                        // prob = probabilities[preferences[v][0]];
                    }

                    double r =
                        ((double)(common_instance->randomInThread(omp_get_thread_num()) %
                                  10000)) /
                        10000.0;
                    if (r <= prob)
                    {
                        if (pref[v] != -1)
                        {
                            return 1;
                        }
                        else
                        {
                            q.push(v);
                            visited[v] = true;
                        }
                    }
                }
            }
        }
    }
    return 0;
}

uint Network::sample_influence_linear_threshold(const kseeds &seeds)
{
    vector<int> pref(number_of_nodes, -1);
    vector<int> c_seeds(Constants::K, 0);

    for (pair<uint, int> p : seeds)
    {
        pref[p.first] = p.second;
        c_seeds[p.second]++;
    }

    int r_node = common_instance->randomInThread(omp_get_thread_num()) % number_of_nodes; // ngấu nhiên 1 node
    if (pref[r_node] != -1)                                                               // if random node is of seeds
    {
        return 1;
    }
    else
    {
        for (int i = 0; i < Constants::K; ++i)
        {
            if (c_seeds[i] > 0)
            {
                queue<uint> q;
                vector<bool> visited(number_of_nodes, false);
                visited[r_node] = true;
                while (true)
                {
                    uint sum_tmp = in_neighbors[r_node].size() * Constants::K;
                    if (sum_tmp == 0)
                        sum_tmp = Constants::K;

                    int r =
                        common_instance->randomInThread(omp_get_thread_num()) % sum_tmp;
                    uint tmp = 0;
                    uint nei = -1;
                    for (uint v : in_neighbors[r_node])
                    {
                        if (tmp <= r && r <= tmp + preferences[v][i] + 1)
                        {
                            nei = v;
                            break;
                        }
                        tmp += (preferences[v][i] + 1);
                    }

                    if (nei == -1)
                        break;
                    else if (pref[nei] == i) // if reach seed node
                    {
                        return 1;
                    }
                    else if (!visited[nei])
                    {
                        r_node = nei;
                        visited[nei] = true;
                    }
                    else
                        break;
                }
            }
        }
    }
    return 0;
}

double Network::my_estimate_influence(const kseeds &seeds)
{
    vector<bool> set(number_of_nodes, false);
    vector<int> set_k(number_of_nodes, 0);
    for (size_t i = 0; i < seeds.size(); ++i)
    {
        set[seeds[i].first] = true;
        set_k[seeds[i].first] = seeds[i].second;
    }
    double val = 0;
    for (int u = 0; u < number_of_nodes; ++u)
    {
        
        vector<pair<int, double>> neis = my_neighbors[u];
        double valU = 0.0;
        for (size_t j = 0; j < neis.size(); ++j)
        {
            int v = neis[j].first;
            if (set[v])
            {
                int k = set_k[v];
                valU += neis[j].second * my_weight_k[v][k];
            }
        }
        valU = pow(valU, alpha[u]);
        val += valU;
    }
    return val;
}

bool Network::read_sensor_data(int no_nodes, string file)
{
    clear();

    this->number_of_nodes = no_nodes;
    ifstream infile(file);

    vector<vector<vector<double>>> tmp(
        no_nodes, vector<vector<double>>(3)); // store all values

    double max_temp = 0.0, min_temp = 1000.0, max_humid = 0.0, min_humid = 1000.0,
           max_light = 0.0, min_light = 10000.0;

    for (string line; getline(infile, line);)
    {
        istringstream iss(line);
        string date, time;
        int epoch, loc;
        double te, hu, li, vol;
        iss >> date >> time >> epoch >> loc >> te >> hu >> li >> vol;
        loc--;

        if (loc >= no_nodes)
            break;

        tmp[loc][0].push_back(te);
        if (te > max_temp)
            max_temp = te;
        if (te < min_temp)
            min_temp = te;

        tmp[loc][1].push_back(hu);
        if (hu > max_humid)
            max_humid = hu;
        if (hu < min_humid)
            min_humid = hu;

        tmp[loc][2].push_back(li);
        if (li > max_light)
            max_light = li;
        if (li < min_light)
            min_light = li;
    }

    vector<double> max_vals = {max_temp, max_humid, max_light};
    vector<double> min_vals = {min_temp, min_humid, min_light};
    vector<int> bin_range = {2, 5, 100};
    vector<int> no_bins(3);
    max_no_bin = 0;
    for (int i = 0; i < 3; ++i)
    {
        no_bins[i] = ceil((max_vals[i] - min_vals[i]) / bin_range[i]);
        if (max_no_bin < no_bins[i])
            max_no_bin = no_bins[i];
    }

    sensor_data = vector<ksensors>(no_nodes, ksensors(3, kbins(max_no_bin, 0.0)));

    for (int i = 0; i < no_nodes; ++i)
    {
        for (int ii = 0; ii < 3; ++ii)
        {
            for (int val : tmp[i][ii])
            {
                int bin_idx = floor((val - min_vals[ii]) / bin_range[ii]);
                sensor_data[i][ii][bin_idx] += (1.0 / tmp[i][ii].size());
            }
        }
    }

    return true;
}

double Network::get_entropy(const kseeds &seeds)
{
    double re = 0.0;

    for (kpoint p : seeds)
    {
        for (double val : sensor_data[p.first][p.second])
            if (val > 0)
                re += (-val * log10(val));
    }
    return re * 100;
}

void Network::clear()
{
    out_neighbors.clear();
    in_neighbors.clear();
    preferences.clear();
    probabilities.clear();
    map_node_id.clear();
    sensor_data.clear();
}

void Network::recursive_entropy(int idx, const kseeds &seeds, double &re, double &prob)
{
    if (idx < seeds.size())
    {
        kpoint p = seeds[idx];
        for (double val : sensor_data[p.first][p.second])
        {
            if (val > 0)
            {
                double prob_new = prob * val;
                recursive_entropy(idx + 1, seeds, re, prob_new);
            }
        }
    }
    else
        re += -prob * log(prob);
}
