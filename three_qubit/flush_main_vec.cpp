#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
//#include <Eigen>
//#include <KroneckerProduct>

#include "omp.h"
#include "basic_funcs_vec.h"
#include "time.h"

using namespace std;
using namespace Eigen;

void printRunTime(time_t t0, time_t t1) {
    int t_run, d, h, m, s; 
    t_run = difftime(t1, t0);
    d = floor(t_run/(24*3600));
    h = floor(t_run/3600);
    m = floor(t_run/60); m = m%60;
    s = t_run%60;
    cout << "TIME TO RUN:\n"
         << d << "d " 
         << h << "h " 
         << m << "m " 
         << s << "s " << endl;
    return;
}

void outputPlotData(string filename, ArrayXXf& FdataList) {
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
/*
void getMatrix(basic_funcs& bf, int Ncycles, int mult, int l, float k, ArrayXf cx, ArrayXf cy, MatrixXcd *rhoList, ArrayXXf FdataList, MatrixXf& matrixM) {
    float collapseOn, collapseOff;
    collapseOn = 1e-3/(k*10); collapseOff = 0.03;

    // basic_funcs bf(collapseOn, collapseOff);
    
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
*/
int main() {
    string evolve_file, pulse_file;
    time_t t0, t1;
    int tp, tf, num_vec, numFidelities, Ncycles, Ohm, listLength, maxIt;
    float dt, dc, acc, collapseOn, collapseOff, J, F;
    bool flush, checking_min;
    Vector2cd s0, s1;
    s0 << 1, 0; s1 << 0, 1;
    VectorXcd ket000, ket111, ket100, ket010, ket001, ket110, ket101, ket011, k000100, finalState;
    VectorXcd ls000[] = {s0,s0,s0,s0,s0,s0}; VectorXcd ls111[] = {s1,s1,s1,s0,s0,s0}; 
    VectorXcd ls100[] = {s1,s0,s0,s0,s0,s0}; VectorXcd ls010[] = {s0,s1,s0,s0,s0,s0};
    VectorXcd ls001[] = {s0,s0,s1,s0,s0,s0}; VectorXcd ls110[] = {s1,s1,s0,s0,s0,s0};
    VectorXcd ls101[] = {s1,s0,s1,s0,s0,s0}; VectorXcd ls011[] = {s0,s1,s1,s0,s0,s0};
    VectorXcd ls000100[] = {s0,s0,s0,s1,s0,s0};
    num_vec = 6;

    ArrayXf cx(20), cy(20); 
    ArrayXf pulse_c[2];
    cx.setZero(); cy.setZero();
    pulse_c[0].setZero(20); pulse_c[1].setZero(20);// pulse_c[2].setZero(20);
    tp = 40; tf = 40; Ncycles = 3; numFidelities = 2;
    dt = 0.01; dc = 0.0001; acc = 1e-5;
    collapseOn = 1e-3/(2*10); collapseOff = 0.03; J = 0.02;
    flush = 1; checking_min = 0;

    ArrayXf fidelities(numFidelities + 1);
    fidelities.setZero();

    basic_funcs_vec bf(collapseOn, collapseOff, J);

    ket000 = bf.tensor_vec(ls000, num_vec); ket111 = bf.tensor_vec(ls111, num_vec);
    ket100 = bf.tensor_vec(ls100, num_vec); ket010 = bf.tensor_vec(ls010, num_vec);
    ket001 = bf.tensor_vec(ls001, num_vec); ket110 = bf.tensor_vec(ls110, num_vec);
    ket101 = bf.tensor_vec(ls101, num_vec); ket011 = bf.tensor_vec(ls011, num_vec);
    k000100 = bf.tensor_vec(ls000100, num_vec);

    // MatrixXcd rhoList[8] = {rho000, rho100, rho010, rho001, rho110, rho101, rho011, rho111};

    // evolve_file = "./outFiles/output_" + to_string(tp) + "_" + to_string(tf) + ".dat";

    maxIt = 0;
    // evolve_file = "./output_vec_" + to_string(maxIt) + ".dat";
    evolve_file = "./output_vec_tp=" + to_string(tp) + ".dat";
    evolve_file = "./output_vec_tp=" + to_string(tp) + ".pls";
    ArrayXXf dataList;
    // if(flush) evolve_file = "./outFiles/yes_coupling.dat";
    // else evolve_file = "./outFiles/no_coupling.dat";

    /*
    evolve_file = "./outFiles/outputF" + to_string(numFidelities) + "_" + to_string(tp);
    if(checking_min) evolve_file += "_min";
    pulse_file = evolve_file + ".pls";
    evolve_file += ".dat";
    */

    // cx << 0.0200494,4.03523e-05,-0.000354767,1.01328e-05,-0.000664413,2.54512e-05,-0.0010761,-0.000144541,-0.000993192,-1.40071e-05,0.000656784,8.40425e-06,2.15173e-05,0.00011009,0.000214219,3.09944e-06,0.000324488,0.000138164,0.000117004,0.000171006;
    // cy << 7.86781e-06,-7.40886e-05,-5.45979e-05,-7.19428e-05,1.00732e-05,0.000286579,-7.75456e-05,0.000314772,-0.000200331,-0.000244141,-2.58088e-05,-0.00012368,-6.00219e-05,-0.00010401,2.69413e-05,-2.68817e-05,-9.77516e-06,0.000133336,-0.00010711,0.00128168;
    // cx[0] = 0.061685; cy[0] = 0.05;
    cx[0] = 0.01;
    pulse_c[0] = cx; pulse_c[1] = cy;
    

    time(&t0);

    cout << "**************** ENTERING TIMED SECTION *****************" << endl;
    int t_cyc[] = {tp, tf};
    Ohm = 0;

    bf.optimizePulse(tp, dt, maxIt, dc, acc, cx, cy, ket000, ket100, ket010, k000100, fidelities, dataList, numFidelities, checking_min);

    // bf.getFidelity(cx, cy, tp, dt, rho000, rho100, rho010, r100100, r100NNN, r010NNN, numFidelities, fidelities, dataList, checking_min);
    
    // bf.evolveState(dt, Ncycles, rho111, rho111, tp, tf, pulse_c, Ohm, flush, dataList, F, finalState);
    outputPlotData(evolve_file, dataList);
    // cout << finalState << endl;

    // ArrayXXf prob00to10[5];
    // ArrayXf optVal;
    // optVal.setZero(5);    
    // for(int k = 1; k < 6; k++) {
    //     optimizeFlushCycle(10, k, 1, 13, rhoList, cx, cy, prob00to10[k - 1]);
    //     optVal[k - 1] = prob00to10[k - 1].rowwise().maxCoeff()(1);
    // }
    
    time(&t1);

    // for(int k = 0; k < 5; k++) cout << optVal[k] << ":\n" << prob00to10[k] << endl << endl;
    printRunTime(t0, t1);
    
    return 0;
}
