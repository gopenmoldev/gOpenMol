/*Copyright 2003 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
#ifndef complex_defined
#include <math.h>
#include <stdio.h>

namespace Plugin {
namespace Fchk {
    
class complex
{
public:
    double real;
    double imaginary;
    complex operator*(complex);
    complex operator^(complex);
    complex operator*(double);
    complex operator/(complex);
    complex operator/(double);
    complex operator+(complex);
    complex operator+(double);
    complex operator-(complex);
    double Re();
    double Im();
    complex Cj();
    void operator*=(complex);
    void operator*=(double);
    void operator^=(complex);
    complex operator/=(complex);
    void operator/=(double);
    void operator+=(complex);
    void operator+=(double);
    void operator-=(complex);
    complex operator=(complex);
    complex operator=(double);
    int operator==(complex);
    int operator==(double);
    int operator!=(complex);
    int operator!=(double);
    complex operator-(void);
    void Argand(double*);
    double magnitude();
    double phase();
    double dphase();
};

complex Complex(double,double);
complex Exp(complex);
double fabs(complex);
 
} // namespace Fchk
} // namespace Plugin

#define complex_defined 1
#endif
