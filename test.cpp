#include <iostream>
#include <cmath>

class Complex {
private:
    double real;
    double imag;

public:
    // Constructors
    Complex(double real = 0.0, double imag = 0.0) : real(real), imag(imag) {}

    // Accessors
    double getReal() const { return real; }
    double getImag() const { return imag; }

    // Addition
    Complex operator+(const Complex& other) const {
        return Complex(real + other.real, imag + other.imag);
    }

    // Subtraction
    Complex operator-(const Complex& other) const {
        return Complex(real - other.real, imag - other.imag);
    }

    // Multiplication
    Complex operator*(const Complex& other) const {
        double real_mult = real * other.real - imag * other.imag;
        double imag_mult = real * other.imag + imag * other.real;
        return Complex(real_mult, imag_mult);
    }

    // Division
    Complex operator/(const Complex& other) const {
        double divisor = other.real * other.real + other.imag * other.imag;
        if (divisor == 0.0) {
            throw std::invalid_argument("Division by zero");
        }
        double real_div = (real * other.real + imag * other.imag) / divisor;
        double imag_div = (imag * other.real - real * other.imag) / divisor;
        return Complex(real_div, imag_div);
    }

    // Exponentiation (using complex pow function)
    Complex operator^(const Complex& exponent) const {
        // Convert to polar form
        double r1 = std::sqrt(real * real + imag * imag);
        double theta1 = std::atan2(imag, real);
        double r2 = std::sqrt(exponent.real * exponent.real + exponent.imag * exponent.imag);
        double theta2 = std::atan2(exponent.imag, exponent.real);

        // Calculate polar power
        double r_pow = std::pow(r1, r2) * std::exp(-r2 * theta1 * std::sin(theta2));
        double theta_pow = r2 * std::sin(theta2) * std::log(r1) + r2 * std::cos(theta2) * theta1;

        // Convert back to rectangular form
        double real_pow = r_pow * std::cos(theta_pow);
        double imag_pow = r_pow * std::sin(theta_pow);
        return Complex(real_pow, imag_pow);
    }

    // Output operator
    friend std::ostream& operator<<(std::ostream& os, const Complex& complex) {
        os << complex.real;
        if (complex.imag >= 0) {
            os << " + " << complex.imag << "i";
        } else {
            os << " - " << -complex.imag << "i";
        }
        return os;
    }
};

int main() {
    // Example usage
    Complex a(1.0, 2.0);
    Complex b(3.0, 4.0);

    std::cout << "Complex a: " << a << std::endl;
    std::cout << "Complex b: " << b << std::endl;

    // Exponentiation
    Complex c = a ^ b;
    std::cout << "a ^ b = " << c << std::endl;

    return 0;
}
