#include "single_qubit.h"

single_qubit::single_qubit() {
//     eye = Matrix2cd::Zero();
//     eyye = Matrix3cd::Zero();
//     I = MatrixXcd::Zero(6, 6);
//     a = Matrix2cd::Zero(); aa = Matrix3cd::Zero();
//     ap = MatrixXcd::Zero(6, 6); apd = MatrixXcd::Zero(6, 6);
//     as = MatrixXcd::Zero(6, 6); asd = MatrixXcd::Zero(6, 6);
//     HX = MatrixXcd::Zero(6, 6); HY = MatrixXcd::Zero(6, 6);
//     HP = MatrixXcd::Zero(6, 6);

//     collapseOn = 0.0; collapseOff = 0.0;
}

single_qubit::single_qubit(int N) {

}

single_qubit::~single_qubit() {

}

float single_qubit::getColOn() const {
    return collapseOn;
}

float single_qubit::getColOff() const {
    return collapseOff;
}


