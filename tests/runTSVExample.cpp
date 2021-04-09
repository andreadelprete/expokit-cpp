#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "utils/readTSV.hpp"

using namespace std;
using namespace Eigen;


int main()
{
    const char fileName[] = "tests/testTSV.tsv";

    // Reading a vector
    vector<MatrixXd> stuff = readTSV(fileName, 9, 1);
    cout << stuff[0] << endl;
    cout << stuff[1] << endl;

    // Reading a matrix
    stuff = readTSV(fileName, 3, 3);
    cout << stuff[0] << endl;
    cout << stuff[1] << endl;

    return 0;
}


