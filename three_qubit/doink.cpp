#include <iostream>
// #include <Eigen>

using namespace std;
// using namespace Eigen;

class parentClass {
    public:
        parentClass();
        parentClass(int, int);
        ~parentClass();
    
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

// ********************* CHILD CLASS******************************

class childClass : public parentClass {
    public:
        childClass();
        childClass(float, int, int);
        ~childClass();

        somefunc(int*, int);
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

childClass::somefunc(int* listOfInts, int sizeOfList) {
    for(int i = 0; i < sizeOfList; i++) cout << ++listOfInts[i] << endl;
}

int main() {
    // parentClass pc(1,2);
    childClass cc(3.4, 7, 5);
    int size = 6;
    int list[size] = {0, 0, 1};

    cc.somefunc(list, size);
    // MatrixXf H;
    // cc.doElse(6, H);

    // cout << H << endl;
    return 0;
}
