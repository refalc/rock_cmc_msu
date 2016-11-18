
#include <iostream>
#include "rock.hpp"
#include <fstream>
#include <ctime>
#include <unistd.h>
#include <sys/time.h>

using namespace std;


int main()
{
    timeval a;
    timeval b;


    dataset set;
    if( ReadSample("SampleTwoDiamonds.txt", set) == -1 ) {
        std::cout << "Error when read sample" << std::endl;
        return -1;
    }

    // todo chose parameters
    double rad = 0.3, threshold = 0.3;
    int number_of_clusters = 2;

    cluster_analysis::rock rock_solver(rad, number_of_clusters, threshold);
    cluster_analysis::cluster_data clusters;

    std::cout << "Start process" << std::endl;
    gettimeofday(&a, 0);
    rock_solver.process(set, clusters);
    gettimeofday(&b, 0);

  
    std::cout << "difference: " << (b.tv_sec - a.tv_sec) << std::endl;
    

    return 0;
}
