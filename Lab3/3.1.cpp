#include <iostream>
#include <math.h>
#include <vector>

double F(double x)
{
    double val = std::sin(x) + x;
    return val;
}

double OmegaVal(double index, std::vector<double>& x)
{
    double val = 1;
    int n = x.size();

    for(int i = 0; i < n; i++)
        if(i == index)
            continue;
        else
            val *= (x[index] - x[i]);
    return val;
}

double NewtonVal(std::vector<double>& P, std::vector<double>& x, double x_star)
{
    int n = P.size();
    double res;
    double value = 0;

    for(int i = 0; i < n; i++)
    {
        res = 1;
        for(int j = 0; j < i; j++)
            res *= (x_star - x[j]);
        value += P[i]*res;
    }
    return value;
}


double LagrangeVal(std::vector<double>& L, std::vector<double>& x, double x_star)
{
    int n = L.size();
    double res;
    double value = 0;

    for(int i = 0; i < n; i++)
    {
        res = 1;
        for(int j = 0; j < n; j++)
            if(j == i)
                continue;
            else
                res *= (x_star - x[j]);
        value += L[i]*res;
    }
    return value;
}

void PrintNewton(std::vector<double>& P, std::vector<double>& x, int n)
{
    for(int i = 0; i < n; i++)
    {
        if(P[i] == 0)
            continue;
        std::cout << P[i]; 
        for(int j = 0; j < i; j++)
        {
            if(x[j] == 0)
                std::cout << "x";
            else if(x[j] > 0)                
                std::cout << "(x - " << x[j] << ")";
            else
               std::cout << "(x + " << -x[j] << ")";
        }
        if(P[i+1] > 0)
            std::cout << " + ";
    }
    std::cout << std::endl;
}

void PrintLagrange(std::vector<double>& L, std::vector<double>& x, int n)
{
    for(int i = 0; i < n; i++)
    {
        if(L[i] == 0)
            continue;
        std::cout << L[i]; 
        for(int j = 0; j < n; j++)
            if(i == j)
                continue;
            else
            {
                if(x[j] == 0)
                    std::cout << "x";
                else if(x[j] > 0)                
                    std::cout << "(x - " << x[j] << ")";
                else
                    std::cout << "(x + " << -x[j] << ")";
            }
            if(L[i+1] > 0)
                std::cout << " + ";
    }
    std::cout << std::endl;
}


double SplitDiff(int i, int j, std::vector<double>& x, std::vector<double>& f)
{
    double val;
    if(j - i == 1)
        val = (f[i] - f[j])/(x[i]-x[j]);
    else
        val = (SplitDiff(i, j-1, x, f) - SplitDiff(i+1, j, x, f))/(x[i]-x[j]);
    return val;
}

int main()
{
    double x_star, val, fval;
    int n;
    //std::cout << "Введите количество интерполяционных узлов: ";
    std::cin >> n;

    std::vector<double> x(n);
    std::vector<double> f(n);
    std::vector<double> L(n);
    std::vector<double> P(n);

    //std::cout << "Введите множество точек - интерполяционных узлов:";
    for(int i = 0; i < n; i++)
    {
        std::cin >> x[i];
        f[i] = F(x[i]);
    }

    //std::cout << "Введите точку, в которой необходимо вычислить значение погрешности интерполяции: ";
    std::cin >> x_star;

    for(int i = 0; i < n; i++)
        L[i] = f[i]/OmegaVal(i, x);

    for(int i = 0; i < n; i++)
        if(i == 0)
            P[i] = f[i];
        else
            P[i] = SplitDiff(0, i, x, f);

    std::cout << "Многочлен Лагранжа:" << std::endl << "L_" << n << "(x) = ";
    PrintLagrange(L, x, n);

    val = LagrangeVal(L, x, x_star);
    fval = F(x_star);
    std::cout << "Значение многочлена Лагранжа в точке X*:" << std::endl << "L_" << n << "(" << x_star << ") = " << val << std::endl;
    std::cout << "Точное значение функции в точке X*:" << std::endl << "f(" << x_star << ") = " << fval << std::endl;
    std::cout << "Значение погрешности интерполяции м-м Лаг-жа: " << std::fabs(fval - val) << std::endl << std::endl;

    std::cout << "Многочлен Ньютона:" << std::endl << "P_" << n << "(x) = ";
    PrintNewton(P, x, n);

    val = NewtonVal(P, x, x_star);
    std::cout << "Значение многочлена Ньютона в точке X*:" << std::endl << "P_" << n << "(" << x_star << ") = " << val << std::endl;
    std::cout << "Точное значение функции в точке X*:" << std::endl << "f(" << x_star << ") = " << fval << std::endl;
    std::cout << "Значение погрешности интерполяции м-м Ньютона: " << std::fabs(fval - val) << std::endl;

    return 0;
}
