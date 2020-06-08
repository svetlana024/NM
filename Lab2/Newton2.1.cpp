#include<iostream>
#include<map>
#include<math.h>
#include<vector>

double F(double x)
{
    double val = pow(3,x) - 5*pow(x,2) + 1; 
    return val;
}

double FDer(double x)
{
    double val = pow(3,x) - 10*x;
    return val;
}

double FDDer(double x)
{
    double val = pow(3,x) - 10;
    return val;
}


bool CheckAB(std::pair<double, double>& p, double eps)
{
    int a, b;
    bool ind = false;
    int fa = F(p.first);
    int fb = F(p.second);

    while(ind == false)
    {
        while(fa*fb > 0)
        {
            if(fa > 0)
            {
                p.first -= eps;
                p.second -= eps;
            }
            else
            {
                p.first += eps;
                p.second += eps;
            }
            fa = F(p.first);
            fb = F(p.second);
        }
            
        if(FDer(p.first)*FDer(p.second) < 0 && FDDer(p.first)*FDDer(p.second) < 0)
        {
            p.first += eps;
            p.second -= eps;
            ind = false;
        }
        else ind = true;
    }
    return true;
}

int main()
{
    double eps;
    int it_count;
    int i = 0;

    std::vector<double> X(2);
    std::vector<std::pair<double, double>> ab(2);
    ab[0] = std::make_pair(0.5, 1.5);
    ab[1] = std::make_pair(3.5, 4.5);

    std::cout << "Введите точность вычислений: ";
    std::cin >> eps;

    for(int i = 0; i < 2; i++)
    {
        CheckAB(ab[i], eps);
        X[1] = ab[i].second;
        it_count = 0;

        while(F(X[1])*FDDer(X[1]) <= 0)
            X[1] -= eps;
        do
        {
            X[0] = X[1];
            X[1] = X[0] - F(X[0])/FDer(X[0]);
            it_count++;
        }
        while(fabs(X[1]-X[0])>= eps);

        std::cout << i + 1 << "-й положительный корень нелинейного уравнения, найденный методом Ньютона: " << X[1] << std::endl;
        std::cout << "Число итераций: " << it_count << std::endl;
    }
    return 0;
}
    
        

    
