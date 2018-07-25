// Operator declarations

#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
// #include <Eigen>
// #include <KroneckerProduct>

using namespace std;
using namespace Eigen;

#define PI 3.14159265358979323846
#define IM std::complex<double> (0.0, 1.0)

#ifndef BASIC_FUNCS_H
#define BASIC_FUNCS_H

class basic_funcs {
    public:
        // Default Constructor
        // basic_funcs();

        // Overload Constructor
        basic_funcs(float, float, float);

        // Destructor
        ~basic_funcs();

        // System functions
        MatrixXcd tensor(MatrixXcd*, int);
        float pulse(float, int, ArrayXf&, int);
        inline void lindbladME(float, float, MatrixXcd&, MatrixXcd&, MatrixXcd&);
        inline void lindbladRK4(float, float, float, MatrixXcd&, MatrixXcd&);
        // void getFidelity(ArrayXf, ArrayXf, int, float, MatrixXcd&, MatrixXcd&, MatrixXcd&, MatrixXcd&, float, float, int, ArrayXf&, ArrayXXf&, int, bool);
        void evolveState(float, int, MatrixXcd&, MatrixXcd&, int*, ArrayXf*, float, bool, ArrayXXf&, float, MatrixXcd&);

        // Member Variables
        float collapseOn, collapseOff;
        
    private:
        // Member Variables
        Matrix2cd a, I, X, Y, Z;
        MatrixXcd X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3,
                  X1S, X2S, X3S, Y1S, Y2S, Y3S, Z1S, Z2S, Z3S,
                  a1, a1d, a2, a2d, a3, a3d,
                  HP, HS, HX1, HX2, HX3, HY1, HY2, HY3;
};

#endif // BASIC_FUNCS_H
