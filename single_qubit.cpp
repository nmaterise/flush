#include "single_qubit.h"

// single_qubit::single_qubit() {
//     eye = Matrix2cd::Zero();
//     eyye = Matrix3cd::Zero();
//     I = MatrixXcd::Zero(6, 6);
//     a = Matrix2cd::Zero(); aa = Matrix3cd::Zero();
//     ap = MatrixXcd::Zero(6, 6); apd = MatrixXcd::Zero(6, 6);
//     as = MatrixXcd::Zero(6, 6); asd = MatrixXcd::Zero(6, 6);
//     HX = MatrixXcd::Zero(6, 6); HY = MatrixXcd::Zero(6, 6);
//     HP = MatrixXcd::Zero(6, 6);

//     collapseOn = 0.0; collapseOff = 0.0;
// }

single_qubit::single_qubit(MatrixXcd tgt, MatrixXcd cst, int NH, int N, float colOn, float colOff) : basic_funcs(NH, N, colOn, colOff) {
    rho00 = kroneckerProduct(ss0, s0);
    rho01 = kroneckerProduct(ss0, s1);
    rho10 = kroneckerProduct(ss1, s0);
    rho11 = kroneckerProduct(ss1, s1);
    rho20 = kroneckerProduct(ss2, s0);
    rho21 = kroneckerProduct(ss2, s1);
    ap = kroneckerProduct(aa, eye); apd = ap.adjoint();
    as = kroneckerProduct(eyye, a); asd = as.adjoint();
    a_ops = {ap, apd, as, asd}
    HX = apd*asd + ap*as; HY = IM*(apd*asd - ap*as);
    HP = -0.1*apd*apd*ap*ap;

    target = tgt; currentState = cst;
}

single_qubit::~single_qubit() {

}

MatrixXcd single_qubit::getTarget() const {
    return target;
}

MatrixXcd single_qubit::getCurrentState() const {
    return currentState;
}

void single_qubit::optimizePulse(float tp, float dt, int maxIt, float dc, float acc, ArrayXf& cx, ArrayXf& cy, MatrixXcd& rho00, MatrixXcd& rho10, MatrixXcd& rho11, MatrixXcd& rho2N, float c1, float c2, int listLength, ArrayXf& fidelities, ArrayXXf& dataList, int numFidelities, bool checking_min) {
    float F, eps;
    int Nmax, it; 
    Nmax = cx.size(); it = 0; eps = dc;
    ArrayXf dFx(Nmax), dFy(Nmax);

    basic_funcs bf(c1, c2);
    
    bf.getFidelity(a_ops, cx, cy, tp, dt, rho00, rho10, rho11, rho2N, c1, c2, fidelities, dataList, listLength, numFidelities, checking_min);
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
            bf.getFidelity(a_ops, cx + eps*diff_vec, cy, tp, dt, rho00, rho10, rho11, rho2N, c1, c2, fidelities, dataList, listLength, numFidelities, checking_min);
            dFx(i) = fidelities(0) - F;
            bf.getFidelity(a_ops, cx, cy + eps*diff_vec, tp, dt, rho00, rho10, rho11, rho2N, c1, c2, fidelities, dataList, listLength, numFidelities, checking_min);
            dFy(i) = fidelities(0) - F;
        }
        cx += dFx; cy += dFy;
        bf.getFidelity(a_ops, cx, cy, tp, dt, rho00, rho10, rho11, rho2N, c1, c2, fidelities, dataList, listLength, numFidelities, checking_min);
        F = fidelities(0);
        it++;
    }
    cout << "=== FINAL FIDELITIES ===\n" << fidelities << endl;
    cout << "=== NUM OF ITERATIONS ===\n" << it << endl;
    return;
}

void single_qubit::getMatrix(int Ncycles, int mult, int l, float k, ArrayXf cx, ArrayXf cy, MatrixXcd *rhoList, ArrayXXf FdataList, MatrixXf& matrixM) {
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

void single_qubit::optimizeFlushCycle(int mult, int k, int mintf, int maxtf, MatrixXcd *rhoList, ArrayXf cx, ArrayXf cy, ArrayXXf& prob00to10) {
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