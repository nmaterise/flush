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
        basic_funcs(int, int, float, float);

        // Destructor
        ~basic_funcs();

        // Accessor Functions
        float getColOn() const;
        float getColOff() const;
        int getHilbertSpace() const;
        int getNumStates() const;

        // System functions
        float pulse(float, int, ArrayXf&, int);
        MatrixXcd Hamiltonian(float, int, ArrayXf& cx, ArrayXf& cy, int, MatrixXcd*);
        // inline void lindbladME(float, float, MatrixXcd&, MatrixXcd&, MatrixXcd&);
        // inline void lindbladRK4(float, float, float, MatrixXcd&, MatrixXcd&, MatrixXcd&);
        inline void lindbladME(float*, MatrixXcd*, MatrixXcd&, MatrixXcd&, MatrixXcd&);
        inline void lindbladRK4(float*, MatrixXcd*, float, MatrixXcd&, MatrixXcd&, MatrixXcd&);
        void getFidelity(MatrixXcd*, ArrayXf, ArrayXf, int, float, MatrixXcd&, MatrixXcd&, MatrixXcd&, MatrixXcd&, float, float, ArrayXf&, ArrayXXf&, int, int, bool);
        void evolveState(float, int, MatrixXcd*, MatrixXcd&, MatrixXcd&, int*, ArrayXf&, ArrayXf&, float, bool, ArrayXXf&, float, MatrixXcd&);

        // Member Variables
        Matrix2cd a, eye, s0, s1;
        Matrix3cd aa, eyye, ss0, ss1, ss2;
        MatrixXcd I;
        
    private:
        // Member Variables
        // Matrix2cd a, s0, s1; 
        // Matrix3cd aa, ss0, ss1, ss2;
        float collapseOn, collapseOff;
        int Hilbert_space, num_states;
};

#endif // BASIC_FUNCS_H
