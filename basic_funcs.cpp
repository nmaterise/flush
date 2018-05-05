#include "basic_funcs.h"

basic_funcs::basic_funcs() {
    eye = Matrix2cd::Zero();
    eyye = Matrix3cd::Zero();
    I = MatrixXcd::Zero(6, 6);
    a = Matrix2cd::Zero(); aa = Matrix3cd::Zero();
    ap = MatrixXcd::Zero(6, 6); apd = MatrixXcd::Zero(6, 6);
    as = MatrixXcd::Zero(6, 6); asd = MatrixXcd::Zero(6, 6);
    HX =MatrixXcd::Zero(6, 6); HY =MatrixXcd::Zero(6, 6);
    HP =MatrixXcd::Zero(6, 6);

    collapseOn = 0.0; collapseOff = 0.0;
}

basic_funcs::basic_funcs(float colOn, float colOff) {
    eye = Matrix2cd::Identity();
    eyye = Matrix3cd::Identity();
    I = MatrixXcd::Identity(6, 6);
    a << 0, 1, 0, 0; aa << 0, 1, 0, 0, 0, sqrt(2), 0, 0, 0;
    ap = kroneckerProduct(aa, eye); apd = ap.adjoint();
    as = kroneckerProduct(eyye, a); asd = as.adjoint();
    HX = apd*asd + ap*as; HY = IM*(apd*asd - ap*as);
    HP = -0.1*apd*apd*ap*ap;

    collapseOn = colOn; collapseOff = colOff;
}

basic_funcs::~basic_funcs() {

}

float basic_funcs::getColOn() const {
    return collapseOn;
}

float basic_funcs::getColOff() const {
    return collapseOff;
}

float basic_funcs::pulse(float t, int tp, ArrayXf& c, int Nmax) {
    float f = 0;
    for(int n = 1; n <= Nmax; n++) {
        f += c[n - 1]*sin(n*PI*t/tp);
    }
    return f;
}

inline void basic_funcs::lindbladME(float col1, float col2, MatrixXcd& rho, MatrixXcd& H, MatrixXcd& output) {
    output = 2.0*IM*PI*(rho*H - H*rho) + col1*(ap*rho*apd - 0.5*(apd*ap*rho + rho*apd*ap)) + col2*(as*rho*asd - 0.5*(asd*as*rho + rho*asd*as));
    return;
}

inline void basic_funcs::lindbladRK4(float col1, float col2, float step, MatrixXcd& rho, MatrixXcd& H, MatrixXcd& output) {
    MatrixXcd t1, t2, t3, k1, k2, k3, k4;
    lindbladME(col1, col2, rho, H, k1);
    t1 = rho + 0.5*step*k1; lindbladME(col1, col2, t1, H, k2);
    t2 = rho + 0.5*step*k2; lindbladME(col1, col2, t2, H, k3);
    t3 = rho + step*k3; lindbladME(col1, col2, t3, H, k4);
    output = rho + step*(k1 + 2*k2 + 2*k3 + k4)/6.0;
    return;
}

void basic_funcs::getFidelity(ArrayXf cx, ArrayXf cy, int tp, float dt, MatrixXcd& rho00, MatrixXcd& rho10, MatrixXcd& rho11, MatrixXcd& rho2N, float c1, float c2, ArrayXf& fidelities, ArrayXXf& dataList, int listLength, int numFidelities, bool checking_min) {
    MatrixXcd currentState1, currentState2, H;
    float F;//, F1, F2, F3;
    int Nmax, Findex;
    Nmax = cx.size();
    currentState1 = rho00; currentState2 = rho10;
    float t = dt;
    for(int i = 0; i < listLength; i++) {
        H = HP + pulse(t, tp, cx, Nmax)*HX + pulse(t, tp, cy, Nmax)*HY;
        lindbladRK4(collapseOn, collapseOff, dt, currentState1, H, currentState1);
        lindbladRK4(collapseOn, collapseOff, dt, currentState2, H, currentState2);
        dataList(0, i) = t;
        dataList(1, i) = (rho11*currentState1.adjoint()).cwiseAbs().trace();
        dataList(2, i) = (rho10*currentState2.adjoint()).cwiseAbs().trace();
        dataList(3, i) = (rho2N*currentState2.adjoint()).cwiseAbs().trace();
        t += dt;
    }
    fidelities(1) = dataList(1, listLength - 1);
    F = fidelities(1);

    if(checking_min) {
        fidelities(2) = dataList(2, listLength - 1);
        fidelities(3) = dataList.rowwise().minCoeff()(2);
    }
    else {
        Findex = floor((listLength - 1)/(numFidelities - 1));
        for(int f = 1; f < numFidelities; f++) {
            fidelities(f + 1) = dataList(2, f*Findex);
            F *= fidelities(f + 1);
        }     
    }
    fidelities(0) = F;

    // F1 = dataList(1, listLength - 1);
    // F2 = dataList(2, listLength - 1);
    // F3 = dataList.rowwise().minCoeff()(2);
    // F3 = dataList(2, floor((listLength - 1)/2));
    // F = F1*F2*F3;
    // fidelities(0) = F;
    // fidelities(1) = F1;
    // fidelities(2) = F2;
    // fidelities(3) = F3;
    return;
}

void basic_funcs::evolveState(float dt, int Ncycles, MatrixXcd& initial, MatrixXcd& target, int *t_cyc, ArrayXf& cx, ArrayXf& cy, float Ohm, bool flush, ArrayXXf& dataList, float F, MatrixXcd& finalState) {
    float collapse;
    int tp, tf, tmax, Nmax, dataIndex, tcycle, tcurrent;
    MatrixXcd currentState, H;
    Nmax = cx.size();
    ArrayXf c1(Nmax), c2(Nmax);
    tp = t_cyc[0]; tf = t_cyc[1];
    currentState = initial; tcurrent = 0; dataIndex = 0;
    finalState = MatrixXcd::Zero(6,6);
    // cout << finalState << endl << "===============" << endl;
    dataList(0, 0) = 0;
    dataList(1, 0) = (target*currentState.adjoint()).cwiseAbs().trace();
    for(int i = 0; i < Ncycles; i++) {
        if(flush) {
            if(i%2 == 0) {
                tcycle = tp; collapse = collapseOn; c1 = cx; c2 = cy;
            } else {
                tcycle = tf; collapse = collapseOff; c1.setZero(); c2.setZero();
            }
        } else {
            tcycle = tp; collapse = collapseOn; c1 = cx; c2 = cy;
        }
        tmax = ceil(tcycle/dt);
        for(int t = 1; t <= tmax; t++) {
            if(flush) H = HP + pulse(t*dt, tcycle, c1, Nmax)*HX + pulse(t*dt, tcycle, c2, Nmax)*HY;
            else H = HP + Ohm*HX;
            dataList(0, dataIndex) = tcurrent + t*dt;
            dataList(1, dataIndex) = (target*currentState.adjoint()).cwiseAbs().trace();
            lindbladRK4(collapseOn, collapse, dt, currentState, H, currentState);
            dataIndex++;
        }
        tcurrent += tcycle;
    }
    int quarterPoint = dataList.cols()*0.25;
    F = dataList.block(0, 3*quarterPoint, 2, quarterPoint).rowwise().mean()(1);
    finalState = currentState;
    // cout << finalState << endl << "===============" << endl;
    return;
}

