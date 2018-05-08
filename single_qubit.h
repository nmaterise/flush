// Operator declarations

#include <iostream>
#include <Eigen>
#include <KroneckerProduct>
#include "basic_funcs.h"

using namespace std;
using namespace Eigen;

#define PI 3.14159265358979323846
#define IM std::complex<double> (0.0, 1.0)

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
        
    private:
        // Member Variables
        Matrix2cd eye, a; 
        Matrix3cd eyye, aa; 
        MatrixXcd I, ap, apd, as, asd, HX, HY, HP;
        float collapseOn, collapseOff;
};

#endif
