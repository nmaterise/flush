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
        basic_funcs();

        // Overload Constructor
        basic_funcs(float, float);

        // Destructor
        ~basic_funcs();

        // Accessor Functions
        float getColOn() const;
        float getColOff() const;

        // System functions
        float pulse(float, int, ArrayXf&, int);
        inline void lindbladME(float, float, MatrixXcd&, MatrixXcd&, MatrixXcd&);
        inline void lindbladRK4(float, float, float, MatrixXcd&, MatrixXcd&, MatrixXcd&);
        void getFidelity(ArrayXf, ArrayXf, int, float, MatrixXcd&, MatrixXcd&, MatrixXcd&, MatrixXcd&, float, float, ArrayXf&, ArrayXXf&, int, int, bool);
        void evolveState(float, int, MatrixXcd&, MatrixXcd&, int*, ArrayXf&, ArrayXf&, float, bool, ArrayXXf&, float, MatrixXcd&);

        // Member Variables
        Matrix2cd eye;
        Matrix3cd eyye;
        MatrixXcd I;
        
    private:
        // Member Variables
        Matrix2cd a; 
        Matrix3cd aa; 
        MatrixXcd ap, apd, as, asd, HX, HY, HP;
        float collapseOn, collapseOff;
};

#endif // BASIC_FUNCS_H
