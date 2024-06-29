#include<iostream>
#include<iomanip>
#include<string>

using namespace std;

class employee
{
    public:
        string name, eid;
        int age;
        char post;

    private:
        float base_salary, increament;
        float current_salary;

    public:
        void get_public_info();
        void get_private_info();
        void make_salary();
};

void employee :: get_public_info()
{
    cout<<"Enter name: ";
    cin>>name;
    cout<<"Enter eid: ";
    cin>>eid;
    cout<<"Enter age: ";
    cin>>age;
    cout<<"Enter post: ";
    cin>>post;
    cout<<endl;
}

void employee :: get_private_info()
{
    cout<<"Enter the base salary: ";
    cin>>base_salary;
    cout<<endl;

    switch(post)
    {
        case 'a':
            increament=5;
            break;
        case 'b':
            increament=10;
            break;
        case 'c':
            increament=15;
            break;
        case 'd':
            increament=20;
            break;
        case 'e':
            increament=25;
            break;
        case 'f':
            increament=30;
            break;
        default:
            increament=0;
            break;

    }
    cout<<endl<<"for post: "<<post<<" increament is: "<<increament<<endl;
}

void employee :: make_salary()
{
    current_salary = base_salary*(1+increament/100);
    cout<<"final salary of "<<name<<" is: "<<current_salary<<endl;
}

int main()
{
    employee sagar;
    sagar.get_public_info();
    sagar.get_private_info();
    sagar.make_salary();
    return(0);
}