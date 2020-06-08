#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>

double ExactVal(double x)
{
    double val = 2 + x + 2*x*log(fabs(x));
    return val;
}

double func(double x)
{
    return 0;    
}

double q_func(double x)
{
    double val = 1/(x*(x-1));
    return val;
}

double p_func(double x)
{
    double val = -1/(x-1);
    return val;
}

void FDM(std::vector<double>& y, std::vector<double>& eps, double h, int l, int r)
{
    int n = (r-l)/h + 1;
    double y_l = 3 + 2*log(2);
    double y_r = -4;
    double nu_0 = 0;
    double nu_1 = 1;
    double mu_0 = 1;
    double mu_1 = -3;
    double a, b, c, d, p, q, f;

    std::vector<double> P(n);
    std::vector<double> Q(n);

    //std::cout << "Вводите ненулевые коэффициенты матрицы и коэффициенты вектора правых частей по порядку: " << std::endl;

    for( int i = 0; i < n; i++)             //прямой ход
    {
        p = p_func(l + h*i);
        q = q_func(l + h*i);
        f = func(l + h*i);
        if( i == 0)
        {
            b = (nu_0 - mu_0/h)/*-2/(h*(2 - p*h)) + q*h/(2 - p*h) + alpha*/;
            c = mu_0/h/*2/(h*(2 - p*h))*/;
            d = y_l/*y_l + h*f/(2 - p*h)*/;

            P[i] = -c/b;
            Q[i] = d/b;
        }
        else if( i != n-1)
        {
            a = 1/pow(h, 2) - p/(2*h);
            b = -2/pow(h, 2) + q;
            c = 1/pow(h, 2) + p/(2*h);
            d = f;
            P[i] = -c/(b + a * P[i-1]);
            Q[i] = (d - a*Q[i-1])/(b + a*P[i-1]);
        }
        else
        {
            a = -mu_1/h/*-2/(h*(2 + p*h))*/;
            b = nu_1 + mu_1/h/*2/(h*(2 + p*h)) - q*h/(2 + p*h) + beta*/;
            d = y_r/*y_r - h*f/(2 + p*h)*/;
            P[i] = 0;
            Q[i] = (d - a*Q[i-1])/(b + a*P[i-1]);
        }
    }

    for( int i = n-1; i >= 0; i--)              //обратный ход
    {
        if( i == n-1)
            y[i] = Q[i];
        else
            y[i] = P[i]*y[i + 1] + Q[i];
        eps[i] = fabs(ExactVal(l + h*i) - y[i]);
        //std::cout << "  x_" << i << " = " << X[i] << std::endl;
    }
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

    FDM(y, eps, h, a, b);
    FDM(rung, eps2, h_r, a, b);

    std::cout << "Решение задачи конечно-разностным методом задано таблично:" << std::endl;
    for(int i = 0; i < n; i++)
    {
        std::cout << "x_" << i + 1 << " = " << a + i*h << "; y_" << i + 1 << " = " << y[i] << std::endl;
        std::cout << "eps = " << eps[i] << ";  Погрешность по Рунге-Ромбергу = " << fabs(y[i]-rung[i/k])/15 <<  std::endl << std::endl;
    }
    return 0;    
}
