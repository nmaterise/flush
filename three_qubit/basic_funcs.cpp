#include "basic_funcs.h"

basic_funcs::basic_funcs(float cOn, float cOff, float J) {
    int num_ops = 6;
    I = Matrix2cd::Identity();
    a << 0, 1, 0, 0; X << 0, 1, 1, 0; Z << 1, 0, 0, -1; Y << 0, IM, -IM, 0;
    MatrixXcd X1List[] = {X,I,I,I,I,I}; MatrixXcd X1SList[] = {I,I,I,X,I,I};
    MatrixXcd Y1List[] = {Y,I,I,I,I,I}; MatrixXcd Y1SList[] = {I,I,I,Y,I,I};
    MatrixXcd Z1List[] = {Z,I,I,I,I,I}; MatrixXcd Z1SList[] = {I,I,I,Z,I,I};
    MatrixXcd X2List[] = {I,X,I,I,I,I}; MatrixXcd X2SList[] = {I,I,I,I,X,I};
    MatrixXcd Y2List[] = {I,Y,I,I,I,I}; MatrixXcd Y2SList[] = {I,I,I,I,Y,I};
    MatrixXcd Z2List[] = {I,Z,I,I,I,I}; MatrixXcd Z2SList[] = {I,I,I,I,Z,I};
    MatrixXcd X3List[] = {I,I,X,I,I,I}; MatrixXcd X3SList[] = {I,I,I,I,I,X};
    MatrixXcd Y3List[] = {I,I,Y,I,I,I}; MatrixXcd Y3SList[] = {I,I,I,I,I,Y};
    MatrixXcd Z3List[] = {I,I,Z,I,I,I}; MatrixXcd Z3SList[] = {I,I,I,I,I,Z};
    MatrixXcd a1List[] = {I,I,I,a,I,I};
    MatrixXcd a2List[] = {I,I,I,I,a,I};
    MatrixXcd a3List[] = {I,I,I,I,I,a};

    X1 = tensor(X1List, num_ops); X1S = tensor(X1SList, num_ops);
    Y1 = tensor(Y1List, num_ops); Y1S = tensor(Y1SList, num_ops);
    Z1 = tensor(Z1List, num_ops); Z1S = tensor(Z1SList, num_ops);
    X2 = tensor(X2List, num_ops); X2S = tensor(X2SList, num_ops);
    Y2 = tensor(Y2List, num_ops); Y2S = tensor(Y2SList, num_ops);
    Z2 = tensor(Z2List, num_ops); Z2S = tensor(Z2SList, num_ops);
    X3 = tensor(X3List, num_ops); X3S = tensor(X3SList, num_ops);
    Y3 = tensor(Y3List, num_ops); Y3S = tensor(Y3SList, num_ops);
    Z3 = tensor(Z3List, num_ops); Z3S = tensor(Z3SList, num_ops);
    a1 = tensor(a1List, num_ops); a1d = a1.adjoint();
    a2 = tensor(a2List, num_ops); a2d = a2.adjoint();
    a3 = tensor(a3List, num_ops); a3d = a3.adjoint();

    HP = -J*(Z1*Z2 + Z2*Z3 + Z1*Z3);
    HS = 2*J*(Z1S + Z2S + Z3S);
    HX1 = X1*X1S; HX2 = X2*X2S; HX3 = X3*X3S;
    HY1 = X1*Y1S; HY2 = X2*Y2S; HY3 = X3*Y3S;

    collapseOn = cOn; collapseOff = cOff;
}

basic_funcs::~basic_funcs() {

}

MatrixXcd basic_funcs::tensor(MatrixXcd* matrix_list, int num_matrices) {
    MatrixXcd output = kroneckerProduct(matrix_list[0], matrix_list[1]).eval();
    for(int i = 2; i < num_matrices; i++) output = kroneckerProduct(output, matrix_list[i]).eval();
    return output;
}

float basic_funcs::pulse(float t, int tp, ArrayXf& c, int Nmax) {
    float f = 0;
    for(int n = 1; n <= Nmax; n++) {
        f += c[n - 1]*sin(n*PI*t/tp);
    }
    return f;
}

inline void basic_funcs::lindbladME(float cp, float cs, MatrixXcd& rho, MatrixXcd& H, MatrixXcd& output) {
    output = 2.0*IM*PI*(rho*H - H*rho) 
             + cp*(X1*rho*X1 - rho) + cs*(a1*rho*a1d - 0.5*(a1d*a1*rho + rho*a1d*a1))
             + cp*(X2*rho*X2 - rho) + cs*(a2*rho*a2d - 0.5*(a2d*a2*rho + rho*a2d*a2))
             + cp*(X3*rho*X3 - rho) + cs*(a3*rho*a3d - 0.5*(a3d*a3*rho + rho*a3d*a3));
    return;
}

inline void basic_funcs::lindbladRK4(float col1, float col2, float step, MatrixXcd& H, MatrixXcd& rho) {
    MatrixXcd t1, t2, t3, k1, k2, k3, k4;
    lindbladME(col1, col2, rho, H, k1);
    t1 = rho + 0.5*step*k1; lindbladME(col1, col2, t1, H, k2);
    t2 = rho + 0.5*step*k2; lindbladME(col1, col2, t2, H, k3);
    t3 = rho + step*k3; lindbladME(col1, col2, t3, H, k4);
    rho = rho + step*(k1 + 2*k2 + 2*k3 + k4)/6.0;
    return;
}

void basic_funcs::getFidelity(ArrayXf cx, ArrayXf cy, int tp, float dt, MatrixXcd& rho000, MatrixXcd& rho100, MatrixXcd& rho010, MatrixXcd& r100100, MatrixXcd& r100NNN, MatrixXcd& r010NNN, int numFidelities, ArrayXf& fidelities, ArrayXXf& dataList, bool checking_min) {
    MatrixXcd currentState1, currentState2, currentState3, H;
    float F, t;
    int Nmax, Findex, listLength;
    Nmax = cx.size();
    currentState1 = rho000; currentState2 = rho100; currentState3 = rho010;
    t = dt;
    listLength = tp/dt;
    dataList.setZero(numFidelities + 2, listLength);
    for(int i = 0; i < listLength; i++) {
        H = HP + pulse(t*dt, tp, cx, Nmax)*HX1
               + pulse(t*dt, tp, cy, Nmax)*HY1 + HS;
        lindbladRK4(0, 0, dt, H, currentState1);
        lindbladRK4(0, 0, dt, H, currentState2);
        lindbladRK4(0, 0, dt, H, currentState3);
        dataList(0, i) = t;
        dataList(1, i) = (r100100*currentState1.adjoint()).cwiseAbs().trace();
        dataList(2, i) = (rho100*currentState2.adjoint()).cwiseAbs().trace();
        dataList(3, i) = (rho010*currentState3.adjoint()).cwiseAbs().trace();
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
    return;
}

void basic_funcs::optimizePulse(float tp, float dt, int maxIt, float dc, float acc, ArrayXf& cx, ArrayXf& cy, MatrixXcd& rho000, MatrixXcd& rho100, MatrixXcd& rho010, MatrixXcd& r100100, MatrixXcd& r100NNN, MatrixXcd& r010NNN, ArrayXf& fidelities, ArrayXXf& dataList, int numFidelities, bool checking_min) {
    float F, eps;
    int Nmax, it; 
    Nmax = cx.size(); it = 0; eps = dc;
    ArrayXf dFx(Nmax), dFy(Nmax);

    getFidelity(cx, cy, tp, dt, rho000, rho100, rho010, r100100, r100NNN, r010NNN, numFidelities, fidelities, dataList, checking_min);
    
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
            getFidelity(cx + eps*diff_vec, cy, tp, dt, rho000, rho100, rho010, r100100, r100NNN, r010NNN, numFidelities, fidelities, dataList, checking_min);
            dFx(i) = fidelities(0) - F;
            getFidelity(cx, cy + eps*diff_vec, tp, dt, rho000, rho100, rho010, r100100, r100NNN, r010NNN, numFidelities, fidelities, dataList, checking_min);
            dFy(i) = fidelities(0) - F;
        }
        cx += dFx; cy += dFy;
        getFidelity(cx, cy, tp, dt, rho000, rho100, rho010, r100100, r100NNN, r010NNN, numFidelities, fidelities, dataList, checking_min);
        F = fidelities(0);
        it++;
    }
    cout << "=== FINAL FIDELITIES ===\n" << fidelities << endl;
    cout << "=== NUM OF ITERATIONS ===\n" << it << endl;
    return;
}

void basic_funcs::evolveState(float dt, int Ncycles, MatrixXcd& initial, MatrixXcd& target, int tp, int tf, ArrayXf* pulse_c, float Ohm, bool flush, ArrayXXf& dataList, float F, MatrixXcd& finalState) {
    float collapse, Ohm1, Ohm2, Ohm3;
    int tmax, Nmax, dataIndex, tcycle, tcurrent, listLength;
    MatrixXcd currentState, H;
    Ohm1 = Ohm; Ohm2 = Ohm; Ohm3 = Ohm;
    Nmax = pulse_c[0].size();
    ArrayXf c1(Nmax), c2(Nmax);// c3(Nmax);
    currentState = initial; tcurrent = 0; dataIndex = 0;
    finalState = MatrixXcd::Zero(6,6);
    if(flush) listLength = (tp*ceil(Ncycles/2.0) + tf*floor(Ncycles/2.0))/dt + 1;
    else listLength = (tp*Ncycles)/dt;
    dataList.setZero(2, listLength);
    dataList(0, 0) = 0;
    dataList(1, 0) = (target*currentState.adjoint()).cwiseAbs().trace();
    for(int i = 0; i < Ncycles; i++) {
        if(flush) {
            if(i%2 == 0) {
                tcycle = tp; collapse = collapseOn; c1 = pulse_c[0]; c2 = pulse_c[1];// c3 = pulse_c[2];
            } else {
                tcycle = tf; collapse = collapseOff; c1.setZero(); c2.setZero();// c3.setZero();
            }
        } else {
            tcycle = tp; collapse = collapseOn; c1 = pulse_c[0]; c2 = pulse_c[1];// c3 = pulse_c[2];
        }
        tmax = ceil(tcycle/dt);
        for(int t = 1; t <= tmax; t++) {
            if(flush) H = HP + pulse(t*dt, tcycle, c1, Nmax)*HX1
                             + pulse(t*dt, tcycle, c2, Nmax)*HY1
                             + pulse(t*dt, tcycle, c1, Nmax)*HX2
                             + pulse(t*dt, tcycle, c2, Nmax)*HY2
                             + pulse(t*dt, tcycle, c1, Nmax)*HX3
                             + pulse(t*dt, tcycle, c2, Nmax)*HY3 + HS;
            else H = HP + Ohm1*HX1 + Ohm2*HX2 + Ohm3*HX3 + HS;
            dataList(0, dataIndex) = tcurrent + t*dt;
            dataList(1, dataIndex) = (target*currentState.adjoint()).cwiseAbs().trace();
            lindbladRK4(collapseOn, collapse, dt, H, currentState);
            dataIndex++;
        }
        tcurrent += tcycle;
    }
    int quarterPoint = dataList.cols()*0.25;
    F = dataList.block(0, 3*quarterPoint, 2, quarterPoint).rowwise().mean()(1);
    finalState = currentState;

    return;
}

