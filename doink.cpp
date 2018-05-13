#include <iostream>
#include <Eigen>

using namespace std;
using namespace Eigen;

class parentClass {
    public:
        parentClass();
        // parentClass(int, int);
        ~parentClass();

        MatrixXf Hamiltonian(int);
        void doSomething(float, MatrixXf&);
    
    private:
        // int num1, num2;
};

parentClass::parentClass() {
    cout << "Constructed Parent Class!" << endl;
}

parentClass::~parentClass() {
    cout << "Deconstructed Parent Class :(" << endl;
}

MatrixXf parentClass::Hamiltonian(int N) {
    return MatrixXf::Identity(N, N);
}

void parentClass::doSomething(float input, MatrixXf& H) {
    H = Hamiltonian(5)*input;
}

// ********************* CHILD CLASS******************************

class childClass : public parentClass {
    public:
        childClass();
        ~childClass();

        MatrixXf Hamiltonian(int);
        void doElse(int, MatrixXf&);
};

childClass::childClass() {
    cout << "Constructed child Class!" << endl;
}

childClass::~childClass() {
    cout << "Deconstructed child Class :(" << endl;
}

MatrixXf childClass::Hamiltonian(int N) {
    return MatrixXf::Random(N, N);
}

void childClass::doElse(int k, MatrixXf& H) {
    doSomething(k*0.1, H);
}

int main() {
    // parentClass pc;
    childClass cc;

    MatrixXf H;
    cc.doElse(6, H);

    cout << H << endl;
}
