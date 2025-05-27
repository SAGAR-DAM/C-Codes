#include <iostream>
#include <string>
using namespace std;

struct Employee {
    string name = "Unknown";
    int age = 30;
    float salary = 50000.0f;
    char post = 'a';

    void display() {
        cout << "Name: " << name << ", Age: " << age
             << ", Salary: " << salary << ", Post: " << post << endl;
    }
};

int main() {
    Employee e1{.name = "Sagar", .post = 'b'};  // Only override 2 values

    e1.display();
    return 0;
}
