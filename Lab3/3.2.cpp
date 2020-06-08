#include <iostream>
#include <math.h>
#include <vector>

void FindC(std::vector<double>& c, std::vector<double>& f, std::vector<double>& h, int n)
{
    std::vector<double> P(n);
    std::vector<double> Q(n);

    for( int i = 1; i < n; i++)             //прямой ход
    {
        if( i == 1)
        {
            P[i] = -h[i]/(2*(h[i-1]+h[i]));   //-c/b;
            Q[i] = 3*((f[i+1]-f[i])/h[i] - (f[i]-f[i-1])/h[i-1])/(2*(h[i-1]+h[i]));                         //d/b;
        }
        else if( i != n-1)
        {
            P[i] = -h[i]/( 2*(h[i-1]+h[i]) + h[i-1]*P[i-1]);                                                                         //-c/(b + a * P[i-1]);
            Q[i] = ( 3*((f[i+1]-f[i])/h[i] - (f[i]-f[i-1])/h[i-1]) - h[i-1]*Q[i-1] )/( 2*(h[i-1]+h[i]) + h[i-1]*P[i-1]);           //(d - a*Q[i-1])/(b + a*P[i-1]);
        }
        else
        {
            P[i] = 0;
            Q[i] = ( 3*((f[i+1]-f[i])/h[i] - (f[i]-f[i-1])/h[i-1]) - h[i-1]*Q[i-1] )/( 2*(h[i-1]+h[i]) + h[i-1]*P[i-1]);       //(d - a*Q[i-1])/(b + a*P[i-1]);
        }
    }

    for( int i = n-1; i >= 1; i--)              //обратный ход
    {
        if( i == n-1)
            c[i] = Q[i];
        else
            c[i] = P[i]*c[i + 1] + Q[i];
    }

}

double SplainValue(double x, double a, double b, double c, double d, double x_star)
{
    double value = a + b*(x_star - x) + c*pow(x_star-x, 2) + d*pow(x_star-x, 3);
    return value;
}


void PrintSplain(double x, double a, double b, double c, double d, int n)
{
    double kaf;

    std::cout << "f(x) = " << a;
    for(int j = 0; j < 3; j++)
    {
        if(j == 0)
            kaf = b;
        else if(j == 1)
            kaf = c;
        else kaf = d;

        if(kaf < 0)
        {
            std::cout << " - " << -kaf;
            if(x < 0)
                std::cout << "(x - " << -x << ")^" << j + 1;
            else if(x > 0)
                std::cout << "(x - " << x << ")^" << j + 1;
            else
                std::cout << "x^" << j + 1;
        }
        else if(kaf > 0)
        {
            std::cout << " + " << kaf;
            if(x < 0)
                std::cout << "(x - " << -x << ")^" << j + 1;
            else if(x > 0)
                std::cout << "(x - " << x << ")^" << j + 1;
            else
                std::cout << "x^" << j + 1;
        }
    }
    std::cout << std::endl;
}


void PrintSplains(std::vector<double>& x, std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, int n)
{
    double kaf;
    for(int i = 0; i < n-1; i++)
    {
        std::cout << "Кубический сплайн на интервале [" << x[i] << ", " << x[i+1] << "] : " << "f(x) = " << a[i];
        for(int j = 0; j < 3; j++)
        {
            if(j == 0)
                kaf = b[i];
            else if(j == 1)
                kaf = c[i];
            else kaf = d[i];

            if(kaf < 0)
            {
                std::cout << " - " << -kaf;
                if(x[i] < 0)
                    std::cout << "(x - " << -x[i] << ")^" << j + 1;
                else if(x[i] > 0)
                    std::cout << "(x - " << x[i] << ")^" << j + 1;
                else
                    std::cout << "x^" << j + 1;
            }
            else if(kaf > 0)
            {
                std::cout << " + " << kaf;
                if(x[i] < 0)
                    std::cout << "(x - " << -x[i] << ")^" << j + 1;
                else if(x[i] > 0)
                    std::cout << "(x - " << x[i] << ")^" << j + 1;
                else
                    std::cout << "x^" << j + 1;
            }
        }
        std::cout << std::endl;
    }
}


int main()
{
    double x_star;
    int n;
    //std::cout << "Введите количество интерполяционных узлов: ";
    std::cin >> n;

    std::vector<double> x(n);
    std::vector<double> f(n);
    std::vector<double> h(n-1, 0);
    std::vector<double> a(n-1, 0);
    std::vector<double> b(n-1, 0);
    std::vector<double> c(n-1, 0);
    std::vector<double> d(n-1, 0);

    //std::cout << "Введите множество точек - интерполяционных узлов:";
    for(int i = 0; i < n-1; i++)
    {
        if(i == 0)
            std::cin >> x[i];
        std::cin >> x[i+1];
        h[i] = x[i+1] - x[i];
    }
    //std::cout << "Введите множество значений функции в узлах интерполяции:";
    for(int i = 0; i < n; i++)
        std::cin >> f[i];

    //std::cout << "Введите точку, в которой необходимо вычислить значение погрешности интерполяции: ";
    std::cin >> x_star;
    FindC(c, f, h, n-1);

    for(int i = 0; i < n-1; i++)
    {
        a[i] = f[i];

        if(i != n-2)
        {
            b[i] = (f[i+1]-f[i])/h[i] - h[i]*(c[i+1]+2*c[i])/3;
            d[i] = (c[i+1]-c[i])/(3*h[i]);
        }
        else
        {
            b[i] = (f[i+1]-f[i])/h[i] -h[i]*c[i]*2/3;
            d[i] = -c[i]/(3*h[i]);
        }
    }
    PrintSplains(x, a, b, c, d, n);

    for(int i = 0; i < n; i++)
    {
        if( x[i] <= x_star && x_star <= x[i+1])
        {
            std::cout << "Точка " << x_star << " принадледит отрезку [" << x[i] << ", " << x[i+1] << "], на этом отрезке функция представляется сплайном: " << std::endl;
            PrintSplain(x[i], a[i], b[i], c[i], d[i], n);
            std::cout << "f(" << x_star << ") = " << SplainValue(x[i], a[i], b[i], c[i], d[i], x_star) << std::endl;
            break;
        }
    }
    return 0;
}
