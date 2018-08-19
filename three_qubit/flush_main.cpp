#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
//#include <Eigen>
//#include <KroneckerProduct>

#include "omp.h"
#include "basic_funcs.h"
#include "time.h"

using namespace std;
using namespace Eigen;

void printRunTime(time_t t0, time_t t1) {
    int t_run, h, m, s; 
    t_run = difftime(t1, t0);
    h = floor(t_run/3600);
    m = floor(t_run/60);
    s = t_run%60;
    cout << "TIME TO RUN:\n" 
         << h << " hr " 
         << (m%60) << " min " 
         << s << " sec " << endl;
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
    int tp, tf, num_ops, numFidelities, Ncycles, Ohm, listLength, maxIt;
    float dt, dc, acc, collapseOn, collapseOff, J, F;
    bool flush, checking_min;
    Matrix2cd I, s0, s1;
    I = Matrix2cd::Identity();
    s0 << 1, 0, 0, 0; s1 << 0, 0, 0, 1;
    MatrixXcd rho000, rho111, rho100, rho010, rho001, rho110, rho101, rho011,
              r000NNN, r111NNN, r100NNN, r010NNN, r001NNN, r110NNN, r101NNN, r011NNN, r000100, finalState;
    MatrixXcd ls000[] = {s0,s0,s0,s0,s0,s0}; MatrixXcd ls111[] = {s1,s1,s1,s0,s0,s0}; 
    MatrixXcd ls100[] = {s1,s0,s0,s0,s0,s0}; MatrixXcd ls010[] = {s0,s1,s0,s0,s0,s0};
    MatrixXcd ls001[] = {s0,s0,s1,s0,s0,s0}; MatrixXcd ls110[] = {s1,s1,s0,s0,s0,s0};
    MatrixXcd ls101[] = {s1,s0,s1,s0,s0,s0}; MatrixXcd ls011[] = {s0,s1,s1,s0,s0,s0};
    MatrixXcd N000[] = {s0,s0,s0,I,I,I}; MatrixXcd N111[] = {s1,s1,s1,I,I,I}; 
    MatrixXcd N100[] = {s1,s0,s0,I,I,I}; MatrixXcd N010[] = {s0,s1,s0,I,I,I};
    MatrixXcd N001[] = {s0,s0,s1,I,I,I}; MatrixXcd N110[] = {s1,s1,s0,I,I,I};
    MatrixXcd N101[] = {s1,s0,s1,I,I,I}; MatrixXcd N011[] = {s0,s1,s1,I,I,I};
    MatrixXcd ls000100[] = {s0,s0,s0,s1,s0,s0};
    num_ops = 6;

    ArrayXf cx(20), cy(20); 
    ArrayXf pulse_c[2];
    cx.setZero(); cy.setZero();
    pulse_c[0].setZero(20); pulse_c[1].setZero(20);// pulse_c[2].setZero(20);
    tp = 40; tf = 40; Ncycles = 3; numFidelities = 2;
    dt = 0.1; dc = 0.0001; acc = 1e-5;
    collapseOn = 1e-3/(2*10); collapseOff = 0.03; J = 0.02;
    flush = 1; checking_min = 0;

    ArrayXf fidelities(numFidelities + 1);
    fidelities.setZero();

    basic_funcs bf(collapseOn, collapseOff, J);

    rho000 = bf.tensor(ls000, num_ops); rho111 = bf.tensor(ls111, num_ops);
    rho100 = bf.tensor(ls100, num_ops); rho010 = bf.tensor(ls010, num_ops);
    rho001 = bf.tensor(ls001, num_ops); rho110 = bf.tensor(ls110, num_ops);
    rho101 = bf.tensor(ls101, num_ops); rho011 = bf.tensor(ls011, num_ops);
    r000NNN = bf.tensor(N000, num_ops); r111NNN = bf.tensor(N111, num_ops);
    r100NNN = bf.tensor(N100, num_ops); r010NNN = bf.tensor(N010, num_ops);
    r001NNN = bf.tensor(N001, num_ops); r110NNN = bf.tensor(N110, num_ops);
    r101NNN = bf.tensor(N101, num_ops); r011NNN = bf.tensor(N011, num_ops);
    r000100 = bf.tensor(ls000100, num_ops);

    // MatrixXcd rhoList[8] = {rho000, rho100, rho010, rho001, rho110, rho101, rho011, rho111};

    // evolve_file = "./outFiles/output_" + to_string(tp) + "_" + to_string(tf) + ".dat";

    maxIt = 0;
    evolve_file = "./outFiles/output_" + to_string(maxIt) + ".dat";
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
    cx[0] = 0.01;
    pulse_c[0] = cx; pulse_c[1] = cy;
    

    time(&t0);

    cout << "**************** ENTERING TIMED SECTION *****************" << endl;
    int t_cyc[] = {tp, tf};
    Ohm = 0;

    bf.optimizePulse(tp, dt, maxIt, dc, acc, cx, cy, rho000, rho100, rho010, r000100, r100NNN, r010NNN, fidelities, dataList, numFidelities, checking_min);

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
