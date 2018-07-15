#include <iostream>
#include <Eigen>
#include <KroneckerProduct>

using namespace std;
using namespace Eigen;

#define PI 3.14159265358979323846
#define IM std::complex<double> (0.0, 1.0)

MatrixXcd tensor(MatrixXcd* matrix_list, int num_matrices) {
	MatrixXcd output = kroneckerProduct(matrix_list[0], matrix_list[1]).eval();
	for(int i = 2; i < num_matrices; i++) output = kroneckerProduct(output, matrix_list[i]).eval();
	return output;
}

void lindbladME(float cp, float cs, MatrixXcd& rho, MatrixXcd& output) {
	Matrix2cd a, I, X, Y, Z;
	MatrixXcd X1, X2, X3, Z1, Z2, Z3,
	          a1, a1d, a2, a2d, a3, a3d, H;
	int num_ops = 6;
	a << 0, 1, 0, 0; X << 0, 1, 1, 0; Z << 1, 0, 0, -1; I << 1, 0, 0, 1;
	MatrixXcd X1List[] = {X,I,I,I,I,I};
	MatrixXcd Z1List[] = {Z,I,I,I,I,I};
	MatrixXcd X2List[] = {I,X,I,I,I,I};
	MatrixXcd Z2List[] = {I,Z,I,I,I,I};
	MatrixXcd X3List[] = {I,I,X,I,I,I};
	MatrixXcd Z3List[] = {I,I,Z,I,I,I};
	MatrixXcd a1List[] = {I,I,I,a,I,I};
	MatrixXcd a2List[] = {I,I,I,I,a,I};
	MatrixXcd a3List[] = {I,I,I,I,I,a};

	X1 = tensor(X1List, num_ops); Z1 = tensor(Z1List, num_ops);
	X2 = tensor(X2List, num_ops); Z2 = tensor(Z2List, num_ops);
	X3 = tensor(X3List, num_ops); Z3 = tensor(Z3List, num_ops);
	a1 = tensor(a1List, num_ops); a1d = a1.adjoint();
	a2 = tensor(a2List, num_ops); a2d = a2.adjoint();
	a3 = tensor(a3List, num_ops); a3d = a3.adjoint();
	H = -0.02*(Z1*Z2 + Z2*Z3 + Z1*Z3);

	output = 2.0*IM*PI*H*(rho*H - H*rho)
			 + cp*(X1*rho*X1 - rho) + cs*(a1*rho*a1d - 0.5*(a1d*a1*rho + rho*a1d*a1))
			 + cp*(X2*rho*X2 - rho) + cs*(a2*rho*a2d - 0.5*(a2d*a2*rho + rho*a2d*a2))
			 + cp*(X3*rho*X3 - rho) + cs*(a3*rho*a3d - 0.5*(a3d*a3*rho + rho*a3d*a3));
	return;
}

int main() {
	Matrix2cd s0, s1;
	s0 << 1, 0, 0, 0; s1 << 0, 0, 0, 1;
	MatrixXcd rho111, mat;
	MatrixXcd ls111[] = {s1,s1,s1,s0,s0,s0};
	int num_ops = 6;
	rho111 = tensor(ls111, num_ops);	
	float cs = 1e-3/20;
	lindbladME(cs, 0.03, rho111, mat);
	cout << ((rho111*mat.adjoint()).trace()) << endl;
	return 0;
}
