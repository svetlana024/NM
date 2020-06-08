#include<iostream>
#include<map>
#include<vector>
#include<math.h>

double Phi1(double x1, double x2)
{
    double val = cos(x2) + 3;
    return val;
}

double Phi2(double x1, double x2)
{
    double val = sin(x1) + 3;
    return val;
}


bool CheckPhi(std::vector<std::pair<double, double>>& p, double eps)
{
    double q;
    while(Phi1(p[0].first, p[0].second) < p[0].first || Phi1(p[0].first, p[0].second) > p[0].second || Phi2(p[1].first, p[1].second) < p[1].first || Phi2(p[1].first, p[1].second) > p[1].second)
    {
            p[0].first += eps;
            p[0].second -= eps;   
            p[0].first += eps;
            p[0].second -= eps;
    }
    return true;
}
         

int main()
{
    double eps, q, kaf;
    int it_count = 0;
    std::vector<std::vector<double>> X(2, std::vector<double> (2, 0));
    std::vector<std::pair<double, double>> ab(2);

    ab[0] = std::make_pair(2, 3);
    ab[1] = std::make_pair(3, 4); 

    std::cout << "Введите точность вычислений: ";
    std::cin >> eps;

    /*if(CheckPhi(ab, eps))
    {
        X[0][1] = (ab[0].first + ab[0].second)/2;
        X[1][1] = (ab[1].first + ab[1].second)/2;
    }
    else std::cerr << "Некорректная функция Phi";*/

    q = std::max(fabs(-sin(ab[1].second)), fabs(-sin(ab[1].first)));
    q = std::max(q, fabs(cos(ab[0].second)));
    q = std::max(q, fabs(cos(ab[0].first)));

    if(q >= 1)
        std::cerr << "Некорректная функция Phi";
    else
        kaf = q/(1-q);
    do
    {
        it_count++;
        X[0][0] = X[0][1];
        X[1][0] = X[1][1];

        X[0][1] = Phi1(X[0][0], X[1][0]);
        X[1][1] = Phi2(X[0][0], X[1][0]);

    }
    while(kaf*std::max(fabs(X[0][1] - X[0][0]), fabs(X[1][1] - X[1][0])) >= eps);

    std::cout << "Метод итераций:" << std::endl;
    std::cout << "q = " << q << std::endl;
    std::cout << "Положительный корень системы нелинейных уравнений:" << "(" << X[0][1] << ", " << X[1][1] << ")" << std::endl;
    std::cout << "Число итераций: " << it_count << std::endl; 

}
