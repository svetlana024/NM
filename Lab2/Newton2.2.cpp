#include<iostream>
#include<map>
#include<vector>
#include<math.h>

double F1(double x1, double x2)
{
    double val = x1 - cos(x2) - 3;
    return val;
}

double F2(double x1, double x2)
{
    double val = x2 - sin(x1) - 3;
    return val;
}

int main()
{
    double eps;
    int it_count = 0;
    std::vector<std::vector<double>> X(2, std::vector<double> (2, 0));
    
    X[0][1] = (2 + 3)/2;
    X[1][1] = (3 + 4)/2;

    std::cout << "Введите точность вычислений: ";
    std::cin >> eps;

    do
    {
        it_count++;
        X[0][0] = X[0][1];
        X[1][0] = X[1][1];

        X[0][1] = X[0][0] - (F1(X[0][0], X[1][0]) - F2(X[0][0], X[1][0])*sin(X[1][0]))/(1 + sin(X[1][0])*cos(X[0][0]));
        X[1][1] = X[1][0] - (F2(X[0][0], X[1][0]) + F1(X[0][0], X[1][0])*cos(X[0][0]))/(1 + sin(X[1][0])*cos(X[0][0]));

    }
    while(std::max(fabs(X[0][1] - X[0][0]), fabs(X[1][1] - X[1][0])) >= eps);

    std::cout << "Метод Ньютона:" << std::endl;
    std::cout << "Положительный корень системы нелинейных уравнений:" << "(" << X[0][1] << ", " << X[1][1] << ")" << std::endl;
    std::cout << "Число итераций: " << it_count << std::endl; 

}
