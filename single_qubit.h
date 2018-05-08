// Operator declarations

#include <iostream>
#include <Eigen>
#include <KroneckerProduct>
#include "basic_funcs.h"

using namespace std;
using namespace Eigen;

// #define PI 3.14159265358979323846
// #define IM std::complex<double> (0.0, 1.0)

#ifndef SINGLE_QUBIT_H
#define SINGLE_QUBIT_H

class single_qubit : public basic_funcs {
    public:
        // Default Constructor
        single_qubit();

        // Overload Constructor
        // single_qubit(float, float);

        // Destructor
        ~single_qubit();

        // Accessor Functions

        // System functions
        void optimizePulse(float, float, int, float, float, ArrayXf&, ArrayXf&, MatrixXcd&, MatrixXcd&, MatrixXcd&, MatrixXcd&, float, float, int, ArrayXf&, ArrayXXf&, int, bool);
        void getMatrix(int, int, int, float, ArrayXf, ArrayXf, MatrixXcd*, ArrayXXf, MatrixXf&);
        void optimizeFlushCycle(int, int, int, int, MatrixXcd*, ArrayXf, ArrayXf, ArrayXXf&);
        
    private:
        // Member Variables
        MatrixXcd rho00, rho01, rho10, rho11, rho2N, rho20, rho21, target, currentState;
        Matrix2cd s0, s1; 
        Matrix3cd ss0, ss1, ss2;
        MatrixXcd ap, apd, as, asd, HX, HY, HP;
        float collapseOn, collapseOff;
};

#endif
