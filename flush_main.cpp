#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
// #include <Eigen>
// #include <KroneckerProduct>

#include "omp.h"
#include "basic_funcs.h"
#include "time.h"

using namespace std;
using namespace Eigen;

void mathematica_outputPlotData(string filename, int listLength, ArrayXXf& dataList) {
    ofstream fout;
    fout.open(filename.c_str());
    fout << "{";
    for(int i = 1; i < dataList.rows(); i++) {
        fout << "{{0,0},";
        for(int j = 0; j < listLength; j++) {
            if(j == (listLength - 1)) fout << "{" << dataList(0, j) << "," << dataList(i, j) << "}";
            else fout << "{" << dataList(0, j) << "," << dataList(i, j) << "},";
        }
        if(i == dataList.rows() - 1) fout << "}";
        else fout << "},\n";
    }
    fout << "};\n";
    fout.close();
    return;
}

void mathematica_outputPulseData(string filename, ArrayXf& cx, ArrayXf& cy) {
    ofstream fout;
    fout.open(filename.c_str());
    fout << "{{";
    for(int i = 0; i < cx.size(); i++) {
        if(i == cx.size() - 1) fout << cx[i] << "},{\n";
        else fout << cx[i] << ",";
    }
    for(int i = 0; i < cy.size(); i++) {
        if(i == cy.size() - 1) fout << cy[i] << "}\n";
        else fout << cy[i] << ",";
    }
    fout << "};\n";
    fout.close();
    return;
}

void outputFPlotData(string filename, ArrayXXf& FdataList) {
    ofstream fout;
    fout.open(filename.c_str());
    for(int i = 0; i < FdataList.rows(); i++) {
        for(int j = 0; j < FdataList.cols(); j++) {
            fout << FdataList(i, j) << " ";
        }
        fout << endl;
    }
    fout.close();
    return;
}

void optimizePulse(float tp, float dt, int maxIt, float dc, float acc, ArrayXf& cx, ArrayXf& cy, MatrixXcd& rho00, MatrixXcd& rho10, MatrixXcd& rho11, MatrixXcd& rho2N, float c1, float c2, int listLength, ArrayXf& fidelities, ArrayXXf& dataList, int numFidelities, bool checking_min) {
    float F, eps;
    int Nmax, it; 
    Nmax = cx.size(); it = 0; eps = dc;
    ArrayXf dFx(Nmax), dFy(Nmax);

    basic_funcs bf(c1, c2);
    
    bf.getFidelity(cx, cy, tp, dt, rho00, rho10, rho11, rho2N, c1, c2, fidelities, dataList, listLength, numFidelities, checking_min);
    F = fidelities(0);
    cout << "=== INITIAL FIDELITIES ===\n" << fidelities << endl;
    while(it <  maxIt && (1 - F) > acc) {
        if((1 - F)/acc > 100) eps = dc;
        else if((1 - F)/acc > 10) eps = 0.1*dc;
        else eps = 0.01*dc;
        #pragma omp parallel for
        for(int i = 0; i < Nmax; i++) {
            ArrayXf diff_vec(Nmax);
            diff_vec.setZero();
            for(int j = 0; j < Nmax; j++) if(i == j) diff_vec(j) = 1;
            bf.getFidelity(cx + eps*diff_vec, cy, tp, dt, rho00, rho10, rho11, rho2N, c1, c2, fidelities, dataList, listLength, numFidelities, checking_min);
            dFx(i) = fidelities(0) - F;
            bf.getFidelity(cx, cy + eps*diff_vec, tp, dt, rho00, rho10, rho11, rho2N, c1, c2, fidelities, dataList, listLength, numFidelities, checking_min);
            dFy(i) = fidelities(0) - F;
        }
        cx += dFx; cy += dFy;
        bf.getFidelity(cx, cy, tp, dt, rho00, rho10, rho11, rho2N, c1, c2, fidelities, dataList, listLength, numFidelities, checking_min);
        F = fidelities(0);
        it++;
    }
    cout << "=== FINAL FIDELITIES ===\n" << fidelities << endl;
    cout << "=== NUM OF ITERATIONS ===\n" << it << endl;
    return;
}

void getMatrix(int Ncycles, int mult, int l, float k, ArrayXf cx, ArrayXf cy, MatrixXcd *rhoList, ArrayXXf FdataList, MatrixXf& matrixM) {
    float collapseOn, collapseOff;
    collapseOn = 1e-3/(k*10); collapseOff = 0.03;
    basic_funcs bf(collapseOn, collapseOff);
    int tp, tf, listLength;
    float dt = 0.1;
    tp = 20; tf = mult*l;
    matrixM = MatrixXf::Zero(6, 6);

    if(Ncycles%2 == 0) listLength = Ncycles*(tp + tf)/(2*dt) + 1;
    else listLength = (tp*ceil(Ncycles/2.0) + tf*floor(Ncycles/2.0))/dt + 1;
    MatrixXcd finalState[6];
    FdataList = ArrayXXf::Zero(2, listLength);
    int t_cyc[] = {tp, tf};
    float dummyF = 0;

    for(int j = 0; j <= 5; j++) {
        for(int i = 0; i <= 5; i++) bf.evolveState(dt, Ncycles, rhoList[i], rhoList[j], t_cyc, cx, cy, 0, 1, FdataList, dummyF, finalState[i]);
        for(int i = 0; i <= 5; i++) matrixM(j, i) = finalState[i](j, j).real();
    }
    return;
}

void optimizeFlushCycle(int mult, int k, int mintf, int maxtf, MatrixXcd *rhoList, ArrayXf cx, ArrayXf cy, ArrayXXf& prob00to10) {
    int size = maxtf - mintf;
    ArrayXXf dummyFdata;
    MatrixXf matrixList[size];
    EigenSolver<MatrixXf> esys[size];
    prob00to10.setZero(2, size);    
    ArrayXf evals;
	VectorXf evecs;
    int evalIndex;
    #pragma omp parallel for
    for(int l = 0; l < size; l++) {
        prob00to10(0, l) = mult*(l + mintf);
        getMatrix(2, mult, l + mintf, k, cx, cy, rhoList, dummyFdata, matrixList[l]);
        esys[l].compute(matrixList[l]);
        evals = esys[l].eigenvalues().real();
        evalIndex = 0;
        for(int i = 0; i < evals.size(); i++) if(evals[i] == evals.maxCoeff()) evalIndex = i;
        evecs = esys[l].eigenvectors().col(evalIndex).real();
        prob00to10(1, l) = ((evecs/evecs.sum())(2));
    }
    return;
}

int main() {
    string evolve_file, pulse_file;
    time_t t0, t1;
    MatrixXcd I, rho00, rho01, rho10, rho11, rho2N, rho20, rho21, target, currentState;
    Matrix2cd eye, s0, s1; 
    Matrix3cd eyye, ss0, ss1, ss2; 
    eye = Matrix2cd::Identity(); eyye = Matrix3cd::Identity();
    s0 << 1, 0, 0, 0; ss0 << 1, 0, 0, 0, 0, 0, 0, 0, 0;
    s1 << 0, 0, 0, 1; ss1 << 0, 0, 0, 0, 1, 0, 0, 0, 0;
    ss2 << 0, 0, 0, 0, 0, 0, 0, 0, 1;
    rho00 = kroneckerProduct(ss0, s0);
    rho01 = kroneckerProduct(ss0, s1);
    rho10 = kroneckerProduct(ss1, s0);
    rho11 = kroneckerProduct(ss1, s1);
    rho20 = kroneckerProduct(ss2, s0);
    rho21 = kroneckerProduct(ss2, s1);
    rho2N = kroneckerProduct(ss2, eye);
    MatrixXcd rhoList[6] = {rho00, rho01, rho10, rho11, rho20, rho21};
    target = kroneckerProduct(ss1, eye);
    I = kroneckerProduct(eyye, eye);
    currentState = rho10;
    ArrayXf cx(20), cy(20); cx.setZero(); cy.setZero(); cx[0] = 0.02;
    int listLength, numFidelities, tp, tf;
    float dt, dc, acc, c1, c2;
    bool checking_min;
    tp = 20; tf = 40; numFidelities = 2; checking_min = 0;
    evolve_file = "./outFiles/outputF" + to_string(numFidelities) + "_" + to_string(tp);
    if(checking_min) evolve_file += "_min";
    pulse_file = evolve_file + ".pls";
    evolve_file += ".dat";
    
    dt = 0.1; dc = 0.0001; acc = 1e-5;
    c1, c2 = 0.0;
    listLength = floor(tp/dt);
    ArrayXf fidelities(numFidelities + 1);
    ArrayXXf dataList(numFidelities + 2, listLength), FdataList;
    fidelities.setZero(); dataList.setZero();

    cx << 0.0200494,4.03523e-05,-0.000354767,1.01328e-05,-0.000664413,2.54512e-05,-0.0010761,-0.000144541,-0.000993192,-1.40071e-05,0.000656784,8.40425e-06,2.15173e-05,0.00011009,0.000214219,3.09944e-06,0.000324488,0.000138164,0.000117004,0.000171006;
    cy << 7.86781e-06,-7.40886e-05,-5.45979e-05,-7.19428e-05,1.00732e-05,0.000286579,-7.75456e-05,0.000314772,-0.000200331,-0.000244141,-2.58088e-05,-0.00012368,-6.00219e-05,-0.00010401,2.69413e-05,-2.68817e-05,-9.77516e-06,0.000133336,-0.00010711,0.00128168;

    ArrayXXf prob00to10[5];
    ArrayXf optVal;
    optVal.setZero(5);
    time(&t0);
    
    for(int k = 1; k < 6; k++) {
        optimizeFlushCycle(10, k, 1, 13, rhoList, cx, cy, prob00to10[k - 1]);
        optVal[k - 1] = prob00to10[k - 1].rowwise().maxCoeff()(1);
    }
    
    time(&t1);

    for(int k = 0; k < 5; k++) cout << optVal[k] << ":\n" << prob00to10[k] << endl << endl;
    cout << "TIME TO RUN :: " << difftime(t1, t0) << endl;
    
    return 0;
}

