#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>

double ExactVal(double x)
{
    double val = x - pow(x, 2) + 1;
    return val;
}

double F(double x, double y1, double y2)
{
    double value = 2*(x*y2 - y1)/(pow(x, 2) + 1);
    return value;
}


void Runge_Kutta(std::vector<double>& y, std::vector<double>& eps, double h, int a, int b)
{
    int n = (b-a)/h + 1;
    double y_h = 1;
    int p = 4;
    std::vector<double> K(p);
    std::vector<double> L(p);
    y[0] = 1;
    eps[0] = fabs(ExactVal(a) - y[0]);
    
    for(int i = 0; i < n-1; i++)
    {
        K[0] = h*y_h;
        L[0] = h*F(a + h*i, y[i], y_h);
        K[1] = h*(y_h + L[0]/2);
        L[1] = h*F(a + h*i + h/2, y[i] + K[0]/2, y_h + L[0]/2);
        K[2] = h*(y_h + L[1]/2);
        L[2] = h*F(a + h*i + h/2, y[i] + K[1]/2, y_h + L[1]/2);
        K[3] = h*(y_h + L[2]);
        L[3] = h*F(a + h*i + h, y[i] + K[2], y_h + L[2]);

        y[i+1] = y[i] + (K[0] + 2*K[1] + 2*K[2] + K[3])/6;
        y_h += (L[0] + 2*L[1] + 2*L[2] + L[3])/6;
        eps[i+1] = fabs(ExactVal(a + h*(i+1)) - y[i+1]);
    }
}

int main()
{
    double k = 0.5;
    double h, h_r, x;
    int n, m;
    double y_h;
    int a = 0;
    int b = 1;

    //std::cout << "Введите значение шага h:";
    std::cin >> h;
    h_r = h*k;
    n = (b-a)/h + 1;
    m = (b-a)/h_r + 1;

    std::vector<double> y(n); 
    std::vector<double> eps(n);
    std::vector<double> rung(m);
    std::vector<double> eps2(m);

    Runge_Kutta(y, eps, h, a, b);
    Runge_Kutta(rung, eps2, h_r, a, b);

    std::cout << "Решение задачи методом Рунге-Кутты задано таблично:" << std::endl;
    for(int i = 0; i < n; i++)
    {
        std::cout << "x_" << i + 1<< " = " << a + i*h << "; y_" << i + 1 << " = " << y[i] << std::endl;
        std::cout << "eps = " << eps[i] << ";  Погрешность по Рунге-Ромбергу = " << fabs(y[i]-rung[i/k])/15 <<  std::endl << std::endl;
    }
    return 0;    

}
