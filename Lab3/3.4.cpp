#include <vector>
#include <math.h>
#include <iostream>

int main()
{
    int n;
    double x_star = 0;
    double y_der, y_doub_der = 0;

    //std::cout << "Введите количество интерполяционных узлов: ";

    std::cin >> n;
    std::vector<double> x(n);
    std::vector<double> y(n);

   //std::cout << "Введите множество точек - интерполяционных узлов:";
    for(int i = 0; i < n; i++)
    {
        std::cin >> x[i];
    }

    //std::cout << "Введите множество значений функции в узлах интерполяции:";
    for(int i = 0; i < n; i++)
    {
        std::cin >> y[i];
    }
    //std::cout << "Введите точку, в которой необходимо вычислить производную: ";
    std::cin >> x_star;

    std::cout << "Производная с первым порядком точности:" << std::endl;
    for(int i = 0; i < n-1; i++)
    {
        if(x_star >= x[i] && x_star <= x[i+1])
        {
            y_der = (y[i+1]-y[i])/(x[i+1]-x[i]);
            y_doub_der = 2*( (y[i+2]-y[i+1])/(x[i+2]-x[i+1]) - y_der )/(x[i+2] - x[i]);
            if(x_star != x[i + 1])
                std::cout << "y'(" << x_star << ") = " << y_der << std::endl;
            else
            {
                std::cout << "Левосторонняя производная y'(" << x_star << ") = " << y_der << std::endl;
                std::cout << "Правосторонняя производная y'(" << x_star << ") = " << (y[i+2]-y[i+1])/(x[i+2]-x[i+1]) << std::endl;
            }
            std::cout << "Производная со вторым порядком точности:" << std::endl;
            std::cout << "y'(" << x_star << ") = " << y_der + 0.5*y_doub_der*(2*x_star - x[i] - x[i+1]) << std::endl; 
            std::cout << "Вторая производная:" << std::endl;
            std::cout << "y''(" << x_star << ") = " << "2*( (" << y[i+2];
            if(y[i+1] < 0)
                std::cout << " + " << -y[i+1] << ")/(" << x[i+2];
            else
                std::cout << " - " << y[i+1] << ")/(" << x[i+2];
            if(x[i+1] < 0)
                std::cout << "+" << -x[i+1] << ") - " <<  "(" << y[i+1];
            else std::cout << "-"<< x[i+1] << ") - " <<  "(" << y[i+1];
            if(y[i] < 0)
                std::cout << "+" << -y[i] << ")/(" << x[i+1];
            else std::cout << "-" << y[i] << ")/(" << x[i+1];
 
            if(x[i] < 0)
                std::cout << "+" << -x[i] << ") )/(" << x[i+2];
            else std::cout << "-" << x[i] << ") )/(" << x[i+2];

            if(x[i] < 0)
                std::cout << "+" << -x[i] << ") = " << y_doub_der << std::endl;
            else std::cout << "-" << x[i] << ") = " << y_doub_der << std::endl;             
            break;
        }
    }

    return 0;
}

