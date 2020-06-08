#include <iostream>
#include <vector>
#include <cmath>

typedef struct
{
    int i;
    int j;
    double value;
} Max;

void InitE(std::vector<std::vector<double>>& U, int n)
{
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            if(i == j)
                U[i][j] = 1;
            else
                U[i][j] = 0;
} 

void MultMatrix(std::vector<std::vector<double>>& F, std::vector<std::vector<double>>& S)
{
    int n = F.size();
    int m = S[0].size();
    double a = 0;
    std::vector<double> H(n);

    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < m; j++)
        {
            for(int k = 0; k < n; k++)
            {
                a += F[j][k]*S[k][i];
            }
            H[j] = a;
            a = 0;
        }
        for(int j = 0; j < n; j++)
            S[j][i] = H[j];
    }
}

void MultMatrix(std::vector<std::vector<double>>& F, std::vector<std::vector<double>>& S, std::vector<std::vector<double>>& Res)
{
    int n = F.size();
    double a = 0;

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            for(int k = 0; k < n; k++)
            {
                a += F[j][k]*S[k][i];
            }
            Res[j][i] = a;
            a = 0;
        }
    }
}

void Transporant(std::vector<std::vector<double>>& Matrix)
{
    double h = 0;
    int n = Matrix.size();

    for(int i = 0; i < n; i++)
        for(int j = i + 1; j < n; j++)
        {
            h = Matrix[i][j];
            Matrix[i][j] = Matrix[j][i];
            Matrix[j][i] = h;
        } 
}   

int JacobiRotations(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& U, double eps)
{
    int n = A.size();
    double fi = 0;
    double a = 0;
    double kr = eps + 1;
    int it_c = 0;
    Max m;
    std::vector<std::vector<double>> Res(n, std::vector<double>(n, 0));
    std::vector<std::vector<double>> SV(n, std::vector<double>(n, 0));

    InitE(SV, n);
    m.value = 0;

    while(kr > eps)
    {
        MultMatrix(SV, U);
        SV = U;
        InitE(U, n);
        m.value = 0;
        kr = 0;
        it_c += 1;

        for(int i = 0; i < n; i++)
            for(int j = i + 1; j < n; j++)
                if( fabs(A[i][j]) > m.value)
                {
                    m.value = fabs(A[i][j]);
                    m.i = i;
                    m.j = j;
                }

        fi = atan(2*A[m.i][m.j]/(A[m.i][m.i] - A[m.j][m.j]))/2;

        U[m.i][m.i] = cos(fi);
        U[m.j][m.j] = cos(fi);
        U[m.i][m.j] = -sin(fi);
        U[m.j][m.i] = sin(fi);

        Transporant(U);
        MultMatrix(U, A, Res);

        Transporant(U);
        MultMatrix(Res, U, A);

        for(int i = 0; i < n; i++)
            for(int j = i + 1; j < n; j++)
                kr += pow(A[i][j], 2);
        kr = pow(kr, 0.5);
    }

    U = SV;
    return it_c;
}

int main()
{

    int n;
    double element, eps;

    std::cout << "Введите размерность матрицы: ";
    std::cin >> n;

    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0));
    
    std::cout << "Введите матрицу A: " << std::endl;

    for( int i = 0; i < n; i++)
        for( int j = 0; j < n; j++)
        {
            std::cin >> A[i][j];
            if(i == j)
                U[i][j] = 1;
        } 
    std::cout << "Введите точность: " << std::endl;
    std::cin >> eps;
    
    std::cout << "Число итераций: " << JacobiRotations(A, U, eps) << std::endl;
    std::cout << "Собственные значения:" << std::endl;
    for(int i = 0; i < n; i++)
        std::cout << "lyambda_" << i + 1 << " = " << A[i][i] << std::endl;
    std::cout << std::endl;
    std::cout << "Собственные векторы:" << std::endl;

    for(int i = 0; i < n; i++)
    {
        std::cout << "x_" << i + 1 << " = ";
        for(int j = 0; j < n; j++)
            std::cout << U[j][i] << " ";
        std::cout << std::endl;
    }

    return 0;
}
