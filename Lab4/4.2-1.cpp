#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>

double ExactVal(double x)
{
    double val = 2 + x + 2*x*log(fabs(x));
    return val;
}

double F(double x, double y1, double y2)
{
    double value = (x*y2 - y1)/(x*(x-1));
    return value;
}

void Runge_Kutta(std::vector<double>& y, std::vector<double>& y2, std::vector<double>& eps, double h, double nu, int a, int b)
{
    int n = (b-a)/h + 1;
    double y_h = 1;
    int p = 4;
    std::vector<double> K(p);
    std::vector<double> L(p);
    y[0] = nu;
    y2[0] = 3 + 2*log(2);
    eps[0] = fabs(ExactVal(2) - y[0]);
    
    for(int i = 0; i < n-1; i++)
    {
        K[0] = h*y2[i];
        L[0] = h*F(a + h*i, y[i], y2[i]);
        K[1] = h*(y2[i] + L[0]/2);
        L[1] = h*F(a + h*i + h/2, y[i] + K[0]/2, y2[i] + L[0]/2);
        K[2] = h*(y2[i] + L[1]/2);
        L[2] = h*F(a + h*i + h/2, y[i] + K[1]/2, y2[i] + L[1]/2);
        K[3] = h*(y2[i] + L[2]);
        L[3] = h*F(a + h*i + h, y[i] + K[2], y2[i] + L[2]);

        y[i+1] = y[i] + (K[0] + 2*K[1] + 2*K[2] + K[3])/6;
        y2[i+1] = y2[i] + (L[0] + 2*L[1] + 2*L[2] + L[3])/6;
        eps[i+1] = fabs(ExactVal(a + h*(i+1)) - y[i+1]);
    }
}

void Print(std::vector<double>& y, std::vector<double>& eps)
{
    double h = 0.1;
    int n = (3-2)/h + 1;
    for(int i = 0; i < n; i++)
    {
        std::cout << "x_" << i + 1 << " = " << 2 + i*h << "; y_" << i + 1 << " = " << y[i] << std::endl;
        std::cout << "eps = " << eps[i] << std::endl;
        /*if(i%k == 0)
            std::cout << ";  Погрешность по Рунге-Ромбергу = " << fabs(y[i]-rung[i])/15 <<  std::endl;
        else std::cout << std::endl;*/
    }
}

void Shoot(std::vector<double>& y, std::vector<double>& eps, double h, int a, int b)
{
    int n = (b-a)/h + 1;
    double epsilon = 0.1;

    std::vector<double> y2(n);
    std::vector<double> Phi(3);
    std::vector<double> nu(3);

    nu[0] = 2;
    nu[1] = 2.8;
    Runge_Kutta(y, y2, eps, h, nu[0], a, b);
    Phi[0] = fabs(y2[n-1] - (y[n-1] + 4)/3);

    Runge_Kutta(y, y2, eps, h, nu[1], a, b);
    Phi[1] = fabs(y2[n-1] - (y[n-1] + 4)/3);
    nu[2] = nu[1] - (nu[1] - nu[0])/(Phi[1] - Phi[0])*Phi[1];
   
    {
        nu[2] = nu[1] - (nu[1] - nu[0])/(Phi[1] - Phi[0])*Phi[1];
        Runge_Kutta(y, y2, eps, h, nu[2], a, b);
        Phi[2] = fabs(y2[n-1] - (y[n-1] + 4)/3);

        nu[0] = nu[1];
        nu[1] = nu[2];  

        Phi[0] = Phi[1];
        Phi[1] = Phi[2];

    }while(Phi[2] > epsilon);
}

int main()
{
    double k = 0.5;
    double h, h_r, x;
    int n, m;
    double y_h;
    int a = 2;
    int b = 3;

    //std::cout << "Введите значение шага h:";
    std::cin >> h;
    h_r = h*k;
    n = (b-a)/h + 1;
    m = (b-a)/h_r + 1;

    std::vector<double> y(n); 
    std::vector<double> eps(n);
    std::vector<double> rung(m);
    std::vector<double> eps2(m);

    Shoot(y, eps, h, a, b);
    Shoot(rung, eps2, h_r, a, b);

    std::cout << "Решение задачи методом стрельбы задано таблично:" << std::endl;
    for(int i = 0; i < n; i++)
    {
        std::cout << "x_" << i + 1 << " = " << a + i*h << "; y_" << i + 1 << " = " << y[i] << std::endl;
        std::cout << "eps = " << eps[i] << ";  Погрешность по Рунге-Ромбергу = " << fabs(y[i]-rung[i/k])/15 <<  std::endl << std::endl;
    }
    return 0;    
}
