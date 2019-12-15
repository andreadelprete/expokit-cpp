// include Eigen before f2c to avoid compilation errors!
#include <Eigen/Core>

#ifdef __cplusplus
extern "C" {

#include "dgpadm.h"
#include "dgexpv.h"
#include "dgchbv.h"
#include "clock.h"

// you need to define MAIN__ function because f2c complains if it doesn't find one
int  MAIN__( ) {  return 0; }

}
#endif

#include <stdio.h>
#include <iostream>

#define MAX_PRINT_N 10
#define PRINT_VECTOR(a) if(a.cols()<=MAX_PRINT_N) std::cout<<#a<<"("<<a.rows()<<"x"<<a.cols()<<"): "<<a.transpose().format(CleanFmt)<<std::endl
#define PRINT_MATRIX(a) if(a.cols()<=MAX_PRINT_N) std::cout<<#a<<"("<<a.rows()<<"x"<<a.cols()<<"):\n"<<a.format(CleanFmt)<<std::endl

class ExpoKit
{
public:
    ExpoKit(int matrix_size);

    void resize(int matrix_size);

    void compute_exp();

};
