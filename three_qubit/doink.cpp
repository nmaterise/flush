#include <cstdarg>
#include <iostream>
// #include <Eigen>
// #include <KroneckerProduct>

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

// ******************************* MAIN ***********************************

// MatrixXf tensor(MatrixXf* matrix_list, int num_matrices) {
//     MatrixXf output = kroneckerProduct(matrix_list[0], matrix_list[1]).eval();
//     for(int i = 2; i < num_matrices; i++) output = kroneckerProduct(output, matrix_list[i]).eval();
//     return output;
// }

int anotherOne(int num1, childClass& cc) {
    int listOfInts[] = {num1, num1, num1};
    cc.somefunc(listOfInts, 3);
    return num1;
}

int main() {
    // parentClass pc(1,2);
    childClass cc(3.4, 7, 5);
    int num;

    num = anotherOne(10, cc);
    cout << num << endl;

    // cout << H << endl;
    return 0;
}
