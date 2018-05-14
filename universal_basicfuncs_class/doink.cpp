#include <iostream>
#include <Eigen>

using namespace std;
using namespace Eigen;

class parentClass {
    public:
        parentClass();
        parentClass(int, int);
        ~parentClass();

        MatrixXf Hamiltonian(int);
        void doSomething(float, MatrixXf&);
    
    private:
        // int num1, num2;
};

parentClass::parentClass() {
    cout << "Default constructor for Parent Class!" << endl;
}

parentClass::parentClass(int num1, int num2) {
    cout << "Overload constructor for Parent Class with arguments:\n" << num1 << ", " << num2 << endl;
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
        childClass(float, int, int);
        // childClass() : parentClass() {};
        // childClass(float) : parentClass(2, 4) { };
        ~childClass();

        MatrixXf Hamiltonian(int);
        void doElse(int, MatrixXf&);
};

childClass::childClass() : parentClass() {
    cout << "Default constructor for child Class!" << endl;
}

childClass::childClass(float num1, int num2, int num3) : parentClass(num2, num3) {
    cout << "Overload constructor for child Class with argument:\n" << num1 << endl;
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
    // parentClass pc(1,2);
    childClass cc(3.4, 7, 5);

    // MatrixXf H;
    // cc.doElse(6, H);

    // cout << H << endl;
}
