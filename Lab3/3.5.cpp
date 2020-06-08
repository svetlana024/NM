#include <vector>
#include <math.h>
#include <iostream>


double YValue(double x)
{
    double value = x/(pow(x, 3) + 8);
    //double value = x/pow(3*x + 4, 2);    
    return value;
}

void Methods(std::vector<double>& RectMethod,  std::vector<double>& TrapezMethod, std::vector<double>& Simpson, double x0, double xk, double h, int index)
{
    double y_val;
    int kaf = 4;
    for(double i = x0; i <= xk; i += h)
    {
        if(i != xk)
            RectMethod[index] += YValue((2*i + h)/2);
        y_val = YValue(i); 
        if(i == x0 || i == xk)
        {       
            TrapezMethod[index] += y_val/2;
            Simpson[index] += y_val;
        }
        else
        {
            TrapezMethod[index] += y_val;
            Simpson[index] += y_val*kaf;
            kaf = (kaf == 4) ? 2 : 4;
        }
    }
    RectMethod[index] *= h;
    TrapezMethod[index] *= h;
    Simpson[index] *= h/3;
}


int main()
{
    double x0, xk, h1, h2;
    double ExR, ExT, ExS;
    double ex_val = -0.006294843240543572;
    std::vector<double> RectMethod(2, 0);
    std::vector<double> TrapezMethod(2, 0);
    std::vector<double> Simpson(2, 0);

    //std::cout << "Введите координаты отрезка интегрирования: ";
    std::cin >> x0 >> xk;

    //std::cout << "Введите значения шагов: ";
    std::cin >> h1 >> h2;
    
    Methods(RectMethod, TrapezMethod, Simpson, x0, xk, h1, 0);
    Methods(RectMethod, TrapezMethod, Simpson, x0, xk, h2, 1);

    std::cout << "Определенный интеграл методом прямоугольников c шагом h = " << h1 << " : " << RectMethod[0] << std::endl; 
    std::cout << "Определенный интеграл методом трапеций c шагом h = " << h1 << " : " << TrapezMethod[0] << std::endl; 
    std::cout << "Определенный интеграл методом Симпсона c шагом h = " << h1 << " : " << Simpson[0] << std::endl; 

    std::cout << "Определенный интеграл методом прямоугольников c шагом h = " << h2 << " : " << RectMethod[1] << std::endl; 
    std::cout << "Определенный интеграл методом трапеций c шагом h = " << h2 << " : " << TrapezMethod[1] << std::endl; 
    std::cout << "Определенный интеграл методом Симпсона c шагом h = " << h2 << " : " << Simpson[1] << std::endl; 

    ExR = RectMethod[0] + (RectMethod[0] - RectMethod[1])/(pow(h2/h1, 2) - 1);
    ExT = TrapezMethod[0] + (TrapezMethod[0] - TrapezMethod[1])/(pow(h2/h1, 2) - 1);
    ExS = Simpson[0] + (Simpson[0] - Simpson[1])/(pow(h2/h1, 4) - 1);

    std::cout << "Уточненные значение методом Рунге-Ромберга-Ричардсона: " << ExR << "; " << ExT << "; " << ExS << std::endl;
    std::cout << "Точное значение интеграла: " << ex_val << std::endl;
    std::cout << "Абсолютная погрешность метода прямоугольников: " << fabs(ex_val - ExR) << std::endl;
    std::cout << "Абсолютная погрешность метода трапеций: " << fabs(ex_val - ExT) << std::endl;
    std::cout << "Абсолютная погрешность метода Симпсона: " << fabs(ex_val - ExS) << std::endl;
        
    return 0;
}

    

