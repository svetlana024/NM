#include<iostream>
#include<map>
#include<math.h>
#include<vector>

double Phi1(double x)
{
    double val = pow((pow(3,x) + 1)/5, 0.5);
    return val;
}

double DPhi1(double x)
{
    double val = pow(3,x)*log(3)/(2*pow(5,0.5)*pow( pow(3,x) + 1, 0.5 ));
    return val;
}

double Phi2(double x)
{
    double val = log(5*pow(x,2)-1)/log(3);
    return val;
}

double DPhi2(double x)
{
    double val = 10*x/( (5*pow(x,2) - 1)*log(3));
    return val;
}

int InitPhi(std::pair<double, double>& p, double eps)
{
    int func_id = 0;

    while(func_id == 0)
    {
        if( Phi1(p.first) >= p.first && Phi1(p.first) <= p.second && Phi1(p.second) >= p.first && Phi1(p.second) <= p.second )
            func_id = 1;
        else if( Phi2(p.first) >= p.first && Phi2(p.first) <= p.second && Phi2(p.second) >= p.first && Phi2(p.second) <= p.second )
            func_id = 2;
        else
        {
            p.first += eps;
            p.second -= eps;
        }
    }
    return func_id;
}

int IterMethod(std::vector<double>& X, double q, double eps, int m)
{
    int it_count = 0;
    double kaf = q/(1-q);

    do
    {
        X[0] = X[1];
        if(m == 1)
            X[1] = Phi1(X[0]);
        else
            X[1] = Phi2(X[0]);
        it_count++;
    }
    while( kaf*fabs(X[1]-X[0]) >= eps);
    return it_count;
}

int main()
{
    double eps, q;
    int it_count;
    int i = 0;

    std::vector<double> X(2);
    std::vector<std::pair<double, double>> ab(2);
    ab[0] = std::make_pair(0.5, 1);
    ab[1] = std::make_pair(3.5, 4.5);

    std::cout << "Введите точность вычислений: ";
    std::cin >> eps;

    for(int i = 0; i < 2; i++)
    {
        it_count = 0;
        if( InitPhi(ab[i], eps) == 1)
        {
            q = std::max(fabs(DPhi1(ab[i].first)), fabs(DPhi1(ab[i].second)));

            while(q >= 1)
            {
                ab[i].first += eps;
                ab[i].second -= eps;
                q = std::max(fabs(DPhi1(ab[i].first)), fabs(DPhi1(ab[i].second)));
            }
            X[1] = (ab[i].first + ab[i].second)/2;
            it_count = IterMethod(X, q, eps, 1);
        }
        else
        {
            q = std::max(fabs(DPhi2(ab[i].first)), fabs(DPhi2(ab[i].second)));

            while(q >= 1)
            {
                ab[i].first += eps;
                ab[i].second -= eps;
                q = std::max(fabs(DPhi2(ab[i].first)), fabs(DPhi2(ab[i].second)));
            }
            X[1] = (ab[i].first + ab[i].second)/2;
            it_count = IterMethod(X, q, eps, 2);
        }
        std::cout << "q = " << q << std::endl;
        std::cout << i+1 << "-й положительный корень нелинейного уравнения, найденный м-м итераций: " << X[1] << std::endl;
        std::cout << "Число итераций: " << it_count << std::endl;        
    }
    return 0;
}
