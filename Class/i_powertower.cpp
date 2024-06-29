#include<iostream>
#include<string>
#include<limits>
#include<cmath>

using namespace std;

class complex
{
    private:
        double real,imag;
        double r;
        float theta;
    public: 
        // Constructor
        complex(double r, double i) : real(r), imag(i) {}
        
        void rect_to_polar()
        {
            r = sqrt(real*real+imag*imag);
            theta = atan2(imag, real);
        }

        void polar_to_rect()
        {
            real = r*cos(theta);
            imag = r*sin(theta);
        }

        double getreal()
        {
            return real;
        }

        double getimag()
        {
            return imag;
        }

        double get_r(complex c)
        {
            c.rect_to_polar();
            return r;
        }

        float get_theta(complex c)
        {
            c.rect_to_polar();
            return theta;
        }

        void getcomplex()
        {
            if(imag>=0) {cout<<"The Complex number is: "<<real<<"+"<<imag<<"j"<<endl;}
            else {cout<<"The Complex number is: "<<real<<imag<<"j"<<endl;}
        }

        void getcomplex_polar()
        {
            rect_to_polar();
            cout<<"The complex number in polar form: "<<r<<" exp("<<theta<<"j)"<<endl;
        }

        static complex addcomplex(complex c1, complex c2)
        {
            complex c_added = complex(c1.real+c2.real,c1.imag+c2.imag);
            return c_added;
        }

        static complex subtractcomplex(complex c1, complex c2)
        {
            complex c_added = complex(c1.real-c2.real,c1.imag-c2.imag);
            return c_added;
        }

        static complex multcomplex(complex c1, complex c2)
        {
            c1.rect_to_polar();
            c2.rect_to_polar();
            double r_mult = c1.r*c2.r;
            float theta_mult = (c1.theta + c2.theta);

            if(theta_mult>M_PI) {theta_mult  = -(2*M_PI-theta_mult);}
            if(theta_mult<M_PI) {theta_mult = 2*M_PI+theta_mult;}
            
            double real_mult = r_mult*cos(theta_mult);
            double imag_mult = r_mult*sin(theta_mult);

            complex c_mult = complex(real_mult, imag_mult);
            return c_mult;
        }

        static complex divcomplex(complex c1, complex c2)
        {
            c1.rect_to_polar();
            c2.rect_to_polar();
            double r_div = c1.r/c2.r;
            float theta_div = (c1.theta - c2.theta);

            if(theta_div>M_PI) {theta_div  = -(2*M_PI-theta_div);}
            if(theta_div<M_PI) {theta_div = 2*M_PI+theta_div;}
            
            double real_div = r_div*cos(theta_div);
            double imag_div = r_div*sin(theta_div);

            complex c_div = complex(real_div, imag_div);
            return c_div;
        }

        static complex complexpow(complex c1, complex c2)
        {
            c1.rect_to_polar();
            c2.rect_to_polar();
            double r1 = c1.r, r2 = c2.r;
            float theta1 = c1.theta, theta2 = c2.theta;

            double r_pow = (pow(r1,r2*cos(theta2)))*exp(-r2*theta1*sin(theta2));
            float theta_pow = r2*sin(theta2)*log(r1)+r2*cos(theta2)*theta1;

            double real_pow = r_pow*cos(theta_pow), imag_pow = r_pow*sin(theta_pow);
            complex c_pow = complex(real_pow, imag_pow);
            c_pow.rect_to_polar();
            
            return c_pow;
        }

};

complex repeatpower(complex c, int n)
{
    if(n<0)
    {
        cout<<"Wrong input; n should not be negative!!"<<endl;
    }
    if(n==0) 
    {
        complex c1 = complex(1,0);
        return c1;
    }
    else 
    {
        complex c1 = c.complexpow(c,repeatpower(c,n-1));
        return c1;
    }
};

int main()
{
    int n;
    cout<<"Enter iteration number: ";
    cin>>n;
    cout<<endl;

    complex c=complex(0.589,-0.458);
    
    for(int i=0;i<=n;i++)
    {
        complex c1 = repeatpower(c,i);
        cout<<"powertower["<<"("<<c.getreal()<<","<<c.getimag()<<"),"<<i<<"]:  ";
        c1.getcomplex();
    }

    std::cin.get();
    return(0);
}