// Operator declarations

#include <iostream>
#include <Eigen>
#include <KroneckerProduct>

using namespace std;
using namespace Eigen;

#define PI 3.14159265358979323846
#define IM std::complex<double> (0.0, 1.0)

#ifndef SINGLE_QUBIT_H
#define SINGLE_QUBIT_H

class single_qubit {
    public:
        // Default Constructor
        single_qubit();

        // Overload Constructor
        single_qubit(float, float);

        // Destructor
        ~single_qubit();

        // Accessor Functions
        float getColOn() const;
        float getColOff() const;

        // System functions
        
    private:
        // Member Variables
};

#endif
