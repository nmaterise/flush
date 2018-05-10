#include <iostream>
#include <Eigen>

using namespace std;
using namespace Eigen;

int main() {
    MatrixXf arr[4];
    Matrix2f a; Matrix3f aa;
    a << 0, 1, 0, 0; aa << 0, 1, 0, 0, 0, sqrt(2), 0, 0, 0;
    cout << a << endl << aa << endl;
    arr[0] = a;
    arr[1] = aa;
    arr[2] = a;
    arr[3] = aa;
    for(int i = 0; i < 4; i++) cout << arr[i] << endl;
    // cout << a << endl;
}
