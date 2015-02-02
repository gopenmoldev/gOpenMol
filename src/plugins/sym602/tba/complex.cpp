#include "complex.h"

namespace Plugin {
namespace Symmetry {

//Defines the complex number class.  It needs a little polishing
/*class complex
{
public:
    double real;
    double imaginary;
    complex operator*(complex);  a*b
    complex operator^(complex);  conj(a)*b
    complex operator*(double);   
    complex operator/(complex);
    complex operator/(double);
    complex operator+(complex);
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
    int operator!=(complex);
    void Argand(double*);   returns mag, arg
    double magnitude();
    double phase();         
    double dphase();
};
*/

complex complex::operator*(complex cp)
{
    complex c1;

    c1.real=real*cp.real-imaginary*cp.imaginary;
    c1.imaginary=real*cp.imaginary+imaginary*cp.real;
    return c1;
}

complex complex::operator^(complex cp)
{
    complex c1;

    c1.real=real*cp.real+imaginary*cp.imaginary;
    c1.imaginary=real*cp.imaginary-imaginary*cp.real;
    return c1;
}

complex complex::operator*(double d)
{
    complex c1;

    c1.real=real*d;
    c1.imaginary=imaginary*d;
    return c1;
}

complex complex::operator/(complex cp)
{
    complex c1;
    double r=real;
    double i=imaginary;
    double mag=cp.magnitude();

    c1.real=(r*cp.real+i*cp.imaginary)/(mag*mag);
    c1.imaginary=(i*cp.real-r*cp.imaginary)/(mag*mag);
    return c1;
}

complex complex::operator/(double d)
{
    complex c1;

    c1.real=real/d;
    c1.imaginary=imaginary/d;
    return c1;
}


complex complex::operator+(complex cp)
{
    complex c1;

    c1.real=real+cp.real;
    c1.imaginary=imaginary+cp.imaginary;
    return c1;
}

complex complex::operator+(double d)
{
    complex c1;

    c1.real=real+d;
    c1.imaginary=imaginary;
    return c1;
}

complex complex::operator-(complex cp)
{
    complex c1;

    c1.real=real-cp.real;
    c1.imaginary=imaginary-cp.imaginary;
    return c1;
}

void complex::operator*=(complex cp)
{
    double r=real;
    double i=imaginary;

    real=r*cp.real-i*cp.imaginary;
    imaginary=r*cp.imaginary+i*cp.real;
}

void complex::operator*=(double d)
{
    real*=d;
    imaginary*=d;
}


void complex::operator^=(complex cp)
{
    double r=real;
    double i=imaginary;

    real=r*cp.real+i*cp.imaginary;
    imaginary=r*cp.imaginary-i*cp.real;
}

complex complex::operator/=(complex cp)
{
    complex c1;
    double mag=cp.magnitude();

    c1.real=(real*cp.real+imaginary*cp.imaginary)/(mag*mag);
    c1.imaginary=(imaginary*cp.real-real*cp.imaginary)/(mag*mag);
    real=c1.real;
    imaginary=c1.imaginary;
    return c1;
}


void complex::operator/=(double d)
{
    real/=d;
    imaginary/=d;
}

void complex::operator+=(complex cp)
{
    real=real+cp.real;
    imaginary=imaginary+cp.imaginary;
}

void complex::operator+=(double d)
{
    real+=d;
}

void complex::operator-=(complex cp)
{
    real=real-cp.real;
    imaginary=imaginary-cp.imaginary;
}

complex complex::operator=(complex cp)
{
    complex rc;

    real=rc.real=cp.real;
    imaginary=rc.imaginary=cp.imaginary;
    return rc;
}

complex complex::operator=(double d)
{
    complex rc;

    real=rc.real=d;
    imaginary=rc.imaginary=0.0;
    return rc;
}

int complex::operator==(complex cp)
{
    if(real!=cp.real)
        return 0;
    if(imaginary!=cp.imaginary)
        return 0;
    return 1;
}

int complex::operator==(double d)
{
    if(imaginary)
        return 0;
    if(real!=d)
        return 0;
    return 1;
}

int complex::operator!=(complex cp)
{
    if(real!=cp.real)
        return 1;
    if(imaginary!=cp.imaginary)
        return 1;
    return 0;
}

int complex::operator!=(double d)
{
    if(imaginary)
        return 1;
    if(real!=d)
        return 1;
    return 0;
}

complex complex::operator-()
{
    complex c;

//  imaginary*=-1;
//  real*=-1;
    c.real=-real;
    c.imaginary=-imaginary;
    return c;
}


double complex::Re()
{
    return real;
}

double complex::Im()
{
    return imaginary;
}

complex complex::Cj()
{
    complex c1;

    c1.real=real;
    c1.imaginary=-imaginary;
    return c1;
}

void complex::Argand(double* op)
{
    op[0]=sqrt(real*real+imaginary*imaginary);
    op[1]=atan2(imaginary,real);
}

double complex::magnitude()
{
    return sqrt(real*real+imaginary*imaginary);
}

double complex::phase()
{
    return atan2(imaginary,real);
}

double complex::dphase()
{
    double arg;

    arg=atan2(imaginary,real);
    return arg*180.0/3.14159265359;
}

complex Complex(double r,double i)
{
    complex c1;

    c1.real=r;
    c1.imaginary=i;
    return c1;
}

complex Exp(complex c)
{
    complex rv;

    rv.real=cos(c.imaginary)*exp(c.real);
    rv.imaginary=sin(c.imaginary)*exp(c.real);
    return rv;
}

double fabs(complex c)
{
    return c.magnitude();
}

} // namespace Symmetry
} // namespace Plugin
