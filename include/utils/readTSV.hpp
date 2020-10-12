#ifndef READTSV_H
#define READTSV_H

#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

vector<MatrixXd> readTSV(const std::string &filename, int rows, int cols, char delimiter='\t')
{
    fstream fin(filename, ios::in); // File pointer
    vector<MatrixXd> result;
    string line;
    string splitted;
    MatrixXd* res;

    while (getline(fin, line)) {
        
        stringstream s(line); // used for breaking words
        int i = 0;
        res = new MatrixXd(rows, cols); // new matrix to contain data

        while (getline(s, splitted, delimiter)) {
            res->operator()(i % rows, i / rows) = stod(splitted);; // So ugly I can't even
            i++;
        }
        if(cols != 1) // Only if it is a matrix. Ugly but works.
            res->transposeInPlace();
        result.push_back(*res);
    }

    return result;
}

#endif