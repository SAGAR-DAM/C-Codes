#include <iostream>
#include <string>
using namespace std;


// Struct that will hold the parameters
struct EmployeeParams {
    string name = "Unknown";
    int age = 0;
    float salary = 0.0f;
    char post = 'd';
};


class Employee {
public:
    string name;
    int age;
    float salary;
    char post;

    // Constructor that takes a struct with all parameters
    Employee(struct EmployeeParams params) {
        name = params.name;
        age = params.age;
        salary = params.salary;
        post = params.post;
    }

    void display() {
        cout << "Name: " << name << ", Age: " << age << ", Salary: " << salary << ", Post: " << post << endl;
    }
};


int main() {
    // You can now pass only the parameters you care about
    EmployeeParams params1;
    params1.name = "Alice";  // Specify only the name
    Employee e1(params1);

    EmployeeParams params2;
    params2.age = 28;  // Specify only the age
    Employee e2(params2);

    EmployeeParams params3;
    params3.salary = 60000.0f;  // Specify only the salary
    params3.name = "Sagar";
    Employee e3(params3);

    e1.display();
    e2.display();
    e3.display();

    return 0;
}
