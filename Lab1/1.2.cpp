#include <iostream>
#include <vector>


void PrintX(std::vector<double>& X, int n)
{
    std::cout << "Решение СЛАУ методом прогонки: " << std::endl;
    for(int i = 0; i < n; i++)
        std::cout << "x_" << i + 1 << " = " << X[i] << std::endl;
}

int main()
{
    int n;
    double a, b, c, d = 0;

    //std::cout << "Введите размерность матрицы: ";
    std::cin >> n;

    std::vector<double> P(n);
    std::vector<double> Q(n);
    std::vector<double> X(n, 0);

    //std::cout << "Вводите ненулевые коэффициенты матрицы и коэффициенты вектора правых частей по порядку: " << std::endl;

    for( int i = 0; i < n; i++)             //прямой ход
    {
        if( i == 0)
        {
            std::cin >> b >> c >> d;
            P[i] = -c/b;
            Q[i] = d/b;
        }
        else if( i != n-1)
        {
            std::cin >> a >> b >> c >> d;
            P[i] = -c/(b + a * P[i-1]);
            Q[i] = (d - a*Q[i-1])/(b + a*P[i-1]);
        }
        else
        {
            std::cin >> a >> b >> d;
            P[i] = 0;
            Q[i] = (d - a*Q[i-1])/(b + a*P[i-1]);
        }
    }

    for( int i = n-1; i >= 0; i--)              //обратный ход
    {
        if( i == n-1)
            X[i] = Q[i];
        else
            X[i] = P[i]*X[i + 1] + Q[i];
        //std::cout << "  x_" << i << " = " << X[i] << std::endl;
    }
    
    PrintX(X, n);
    
    return 0;
}
