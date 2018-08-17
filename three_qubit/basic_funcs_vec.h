// Operator declarations

#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
// #include <Eigen>
// #include <KroneckerProduct>

using namespace std;
using namespace Eigen;

#define PI 3.14159265358979323846
#define IM std::complex<double> (0.0, 1.0)

#ifndef BASIC_FUNCS_H_VEC
#define BASIC_FUNCS_H_VEC

class basic_funcs_vec {
    public:
        // Default Constructor
        // basic_funcs_vec();

        // Overload Constructor
        basic_funcs_vec(float, float, float);

        // Destructor
        ~basic_funcs_vec();

        // System functions
        MatrixXcd tensor(MatrixXcd*, int);
        VectorXcd tensor_vec(VectorXcd*, int);
        float pulse(float, int, ArrayXf&, int);
        inline void timePropagator(float, VectorXcd&, MatrixXcd&, VectorXcd&);
        inline void timePropagatorRK4(float, MatrixXcd&, VectorXcd&);
        void getFidelity(ArrayXf, ArrayXf, int, float, VectorXcd&, VectorXcd&, VectorXcd&, VectorXcd&, int, ArrayXf&, ArrayXXf&, bool);
        void optimizePulse(float, float, int, float, float, ArrayXf&, ArrayXf&, VectorXcd&, VectorXcd&, VectorXcd&, VectorXcd&, ArrayXf&, ArrayXXf&, int, bool);
        void evolveState(float, int, MatrixXcd&, MatrixXcd&, int, int, ArrayXf*, float, bool, ArrayXXf&, float, MatrixXcd&);

        // Member Variables
        float collapseOn, collapseOff;
        
    private:
        // Member Variables
        Matrix2cd a, I, X, Y, Z;
        MatrixXcd X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3,
                  X1S, X2S, X3S, Y1S, Y2S, Y3S, Z1S, Z2S, Z3S,
                  a1, a1d, a2, a2d, a3, a3d,
                  HP, HS, HX1, HX2, HX3, HY1, HY2, HY3,
                  Eye;
};

#endif // BASIC_FUNCS_H
