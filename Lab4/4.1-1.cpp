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


void ExEuler(std::vector<double>& y, std::vector<double>& eps, double h, int a, int b)
{
    double n = (b - a)/h + 1;
    double y_h = 1;
    y[0] = 1;
    eps[0] = fabs(ExactVal(a) - y[0]);
    
    for(int i = 0; i < n-1; i++)
    {
        y[i+1] = y[i] + h*y_h;
        y_h += h*F(a + h*i, y[i], y_h);
        eps[i+1] = fabs(ExactVal(a + h*(i+1)) - y[i+1]);
    }
}

int main()
{
    double k = 0.5;
    double h, n, m, x;
    double y_h;
    int b = 1;
    int a = 0;

    //std::cout << "Введите значение шага h:";
    std::cin >> h;
    n = (b-a)/h + 1;
    m = (b-a)/(h*k) + 1;

    std::vector<double> y(n); 
    std::vector<double> eps(n);
    std::vector<double> rung(m);
    std::vector<double> eps2(m);

    ExEuler(y, eps, h, a, b);
    ExEuler(rung, eps2, k*h, a, b);

    std::cout << "Решение задачи явным методом Эйлера задано таблично:" << std::endl;
    for(int i = 0; i < n; i++)
    {
        std::cout << "x = " << a + i*h << "; y = " << y[i] << std::endl;
        std::cout << "eps = " << eps[i] << ";  Погрешность по Рунге-Ромбергу = " << fabs(y[i]-rung[i/k]) <<  std::endl << std::endl;
    }
    return 0;    

}
