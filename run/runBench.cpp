#include <Eigen/Core>
#include <cstdlib>
#include <iostream>

#include "LDS2OrderUtility.hpp"

#include "utils/readTSV.hpp"
#include "utils/stop-watch.h"

using namespace std;
using namespace expokit;

#define MANY 100

#define N 4
#define M N * 3 * 2

/*
    argv[0] - executable name
    argv[1] - minSquarings
*/
int main()
{
    LDS2OrderUtility<double, Dynamic> utilDynamic(M);
    LDS2OrderUtility<double, M> utilStatic();
}