#include <fstream>
#include <iostream>
#include <omp.h>
#include <stdlib.h>
#include <string>
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>
#include <tuple>

#include "Constants.h"
#include "Degree.h"
#include "Greedy.h"
#include "Network.h"
#include "Streaming.h"

//#include <time.h>
#include <chrono>
#include <iomanip> // setprecision
#if defined(_WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

using namespace std;

#pragma warning(disable : 4996)
vector<uint> B_vals = {15, 20, 25, 30, 35, 40, 45, 50};
vector<uint> M_vals = {3};
vector<double> Delta_vals = {1};

Network *g;

void print_help()
{
    cout << "Options: " << endl;
    cout << "-f <data filename>" << endl
         << "-V <size of V> # default: 50" << endl
         << "-t <type of experiment, 0: influence maximization, 1: sensor "
            "placement> # default: 0"
         << endl
         << "-k <value of k> # default: 3" << endl
         << "-B <value of B> # default: 10" << endl
         << "-M <value of M> # default: 3" << endl
         << "-e <value of epsilon> # default: 0.5 for IM and 0 for sensor "
            "placement"
         << endl
         << "-n <value of eta - denoise step for RStream> # default: 2" << endl
         << "-g <value of gamma> # default: 1.0" << endl
         << "-a <algorithm, 1: DStream, 2: RStream, 3: SGr>, 4: SampleRstream # "
            "default: 1, please use SSA source code for testing IM algorithm"
         << endl
         << "-p <number of threads to running algorithms> # default: 4" << endl;
}

pair<string, int> parseArgs(int argc, char **argv)
{
    string filename = "error";
    int no_nodes = 50;

    int i = 1;
    while (i < argc - 1)
    {
        string arg = argv[i];
        if (arg == "-f")
        {
            filename = argv[i + 1];
        }
        else if (arg == "-g" || arg == "-e")
        {
            string s_val = argv[i + 1];
            std::string::size_type sz;
            Constants::EPS = std::stod(s_val, &sz);
        }
        else
        {
            string s_val = argv[i + 1];
            std::string::size_type sz;
            int val = std::stoi(s_val, &sz);
            if (arg == "-t")
            {
                Constants::DATA = val == 0 ? Social : Sensor;
                if (Constants::DATA == Sensor)
                {
                    Constants::NO_DENOISE_STEPS = 1;
                }
            }
            else if (arg == "-V")
            {
                no_nodes = val;
            }
            else if (arg == "-k")
            {
                Constants::K = val;
            }
            else if (arg == "-r")
            {
                Constants::RHO = val;
            }
            else if (arg == "-b")
            {
                Constants::BETA = val;
            }
            else if (arg == "-B")
            {
                Constants::BUDGET = val;
            }
            else if (arg == "-n")
            {
                Constants::NO_DENOISE_STEPS = val;
                if (Constants::DATA == Sensor)
                {
                    Constants::NO_DENOISE_STEPS = 1;
                }
            }
            else if (arg == "-a")
            {
                switch (val)
                {
                case 0:
                    Constants::ALGORITHM = Adaptive;
                    break;
                case 1:
                    Constants::ALGORITHM = eAppro; //Alg2
                    break;
                case 2:
                    Constants::ALGORITHM = ream;
                    break;
                case 3:
                    Constants::ALGORITHM = SGr;
                    break;
                case 4:
                    Constants::ALGORITHM = SampleRstream;
                    break;
                case 5:
                    Constants::ALGORITHM = Tang;
                    break;
                case 6:
                    Constants::ALGORITHM = Tang2;
                    break;
                case 8:
                    Constants::ALGORITHM = eApproPlush;
                    break;
                case 9:
                    Constants::ALGORITHM = Testmpi;
                    break;
                case 10:
                    Constants::ALGORITHM = DS;
                    break;
                case 11:
                    Constants::ALGORITHM = RS;
                    break;
                default:
                    Constants::ALGORITHM = Adaptive;
                    break;
                }
            }
            else if (arg == "-p")
            {
                Constants::NUM_THREAD = val;
            }
        }
        i += 2;
    }

    return pair<string, int>(filename, no_nodes);
}

void run_command(string filename, int no_nodes)
{
    Network *g = new Network();
    kseeds seedsf;
    bool r = false;
    if (Constants::DATA == Social)
    {
        r = g->read_network_from_file2(filename, false);
    }
    else
    {
        r = g->read_sensor_data(no_nodes, filename);
    }

    if (!r)
    {
        cout << "Wrong file, format or arguments" << endl;
        print_help();
    }
    else
    {
        double sol = 0.0;
        int no_queries = 0;
        int size_seedsf = 0;
        omp_set_num_threads(Constants::NUM_THREAD);
        switch (Constants::ALGORITHM)
        {
        case Adaptive:
        {
            Greedy *gr = new Greedy(g);
            estimate e_max(g->get_no_nodes());
            kpoint fe_max;
            double fe;
            sol = gr->get_solution(seedsf, true, e_max, fe_max, fe);
            no_queries = gr->get_no_queries();
            delete gr;
            break;
        }
        case eAppro:
        {
            Greedy *gr = new Greedy(g);
            sol = gr->get_solution2(seedsf, true);
            no_queries = gr->get_no_queries();
            delete gr;
            break;
        }
        case eApproPlush:
        {
            Greedy *gr = new Greedy(g);
            sol = gr->get_solution3(seedsf, true);
            no_queries = gr->get_no_queries();
            delete gr;
            break;
        }
        case RS:
        {
            Constants::NO_DENOISE_STEPS = 1;
            Constants::ALPHA = 2 / (3 + Constants::BETA - (Constants::BETA / Constants::K));
            Constants::EPS_TAG = (3 + Constants::BETA - (Constants::BETA / Constants::K)) * Constants::EPS;
            Streaming *ra = new Streaming(g);
            sol = ra->get_solution(false);
            no_queries = ra->get_no_queries();
            delete ra;
            break;
        }
        case DS:
        {
            Constants::NO_DENOISE_STEPS = 1;
            Constants::ALPHA = 0.5;
            Constants::EPS_TAG = 4 * Constants::EPS;
            Streaming *de = new Streaming(g);
            sol = de->get_solution(true);
            no_queries = de->get_no_queries();
            delete de;
            break;
        }      
        case Greedy:
        {
            Greedy *gr = new Greedy(g);
            sol = gr->get_solutionT2(seedsf);
            no_queries = gr->get_no_queries();
            delete gr;
            break;
        }
        }
        size_seedsf = seedsf.size();

        if (Constants::DATA == Sensor)
            sol = sol / 100;
        switch (Constants::ALGORITHM)
        {
        case Adaptive:
            std::cout << "Linear Apdaptive"
                      << "," << Constants::EPS << ","
                      << Constants::BUDGET << "," << no_nodes << "," << sol << "," << size_seedsf << "," << no_queries << ",";
            break;
        case eAppro:
            std::cout << "eps-approximation"
                      << "," << Constants::EPS << ","
                      << Constants::BUDGET << "," << no_nodes << "," << sol << "," << size_seedsf << "," << no_queries << ",";
            break;
        case eApproPlush:
            std::cout << "eps-approximation++"
                      << "," << Constants::EPS << ","
                      << Constants::BUDGET << "," << no_nodes << "," << sol << "," << size_seedsf << "," << no_queries << ",";
            break;
        case DS:
            std::cout << "DS"
                      << "," << Constants::EPS << ","
                      << Constants::BUDGET << "," << no_nodes << "," << sol << "," << size_seedsf << "," << no_queries << ",";
            break;
        case RS:
            std::cout << "RS"
                      << "," << Constants::EPS << ","
                      << Constants::BUDGET << "," << no_nodes << "," << sol << "," << size_seedsf << "," << no_queries << ",";
            break;       
        case Greedy:
            std::cout << "Greedy"
                      << "," << Constants::EPS << ","
                      << Constants::BUDGET << "," << no_nodes << "," << sol << "," << size_seedsf << "," << no_queries << ",";
            break;
        default:
            std::cout << "Default"
                      << ", epsilon: " << Constants::EPS << ", with Budget: "
                      << Constants::BUDGET << ", f(.)=" << sol << ", num queries:" << no_queries << ", time:";
            break;
        }
    }
    delete g;
}

void print_performances(struct rusage *performance_start, struct rusage *performance_end)
{
    cout << (performance_end->ru_utime.tv_sec -
             performance_start->ru_utime.tv_sec) +
                1e-6 * (performance_end->ru_utime.tv_usec -
                        performance_start->ru_utime.tv_usec)
         << endl;
}

int main(int argc, char *argv[])
{
    srand(time(NULL));
    pair<string, int> r = parseArgs(argc, argv);
    struct rusage performance_start, performance_end;
    getrusage(RUSAGE_SELF, &performance_start);

    if (r.first == "error")
    {
        std::cout << "Wrong file, format or arguments" << endl;
        print_help();
    }
    else
    {
        run_command(r.first, r.second);
    }

    getrusage(RUSAGE_SELF, &performance_end);
    print_performances(&performance_start, &performance_end);
    return 0;
}
