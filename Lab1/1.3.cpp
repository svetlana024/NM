#include <iostream>
#include <vector>
#include <cmath>

void Swap(std::vector<std::vector<double>>& M, int i)
{
    int n = M.size();
    int k = i + 1;

    while(M[k][i] == 0)
        k++;
    std::vector<double> N(n);
    N = M[i];
    M[i] = M[k];
    M[k] = N;
}

void PrintMatrix(std::vector<std::vector<double>>& M)
{
    for(int i = 0; i < M.size(); i++){
        for(int j = 0; j < M[0].size(); j++)
            std::cout << M[i][j] <<  " ";
        std::cout << std::endl;
    }
}

int IterMethod(std::vector<std::vector<double>>& A, std::vector<double>& X, double eps)
{
    int n = A.size();
    int it_count = 0;	
    double norma = 0;
    double kaf = 0;
    double sum = 0;

    std::vector<double> X_1(n, 0);		//вспомогательный вектор, хранит результат пред-й итерации
    
    for( int i = 0; i < n; i++)			// приводим матрицу к виду x = beta + alpha*X
    {
        if( A[i][i] == 0)
            Swap(A, i);
        for( int j = 0; j <= n; j++)
        {
            if(j == i)
                continue;
            else
            {   
                if( j != n)
                {
                    A[i][j] = (-1)* A[i][j]/ A[i][i];
                    sum += fabs(A[i][j]);
                }
                else
                    A[i][j] /= A[i][i];
            }
        }

        if( sum > norma)				//норма матрицы alfha
        {
            norma = sum;
            sum = 0;
        }
        A[i][i] = 0;
    }
    norma /= (1 - norma);				//коэффициент

    for(int i = 0; i < n; i++)
        X_1[i] = A[i][n];

    while(true)
    {
        for(int i = 0; i < n; i++)
        {
            X[i] = 0;
            for(int j = 0; j < n; j++)
                X[i] += A[i][j]* X_1[j];
            X[i] += A[i][n];
            if( kaf < fabs(X[i] - X_1[i]))
		        kaf = fabs(X[i] - X_1[i]);
        }
        it_count++;
        if(eps > kaf*norma)
            break;
        X_1 = X;
        kaf = 0;
    }
    return it_count;
}

int main()
{
    int n;
    double element, eps;

    //std::cout << "Введите размерность матрицы: ";
    std::cin >> n;

    std::vector<std::vector<double>> A(n, std::vector<double>(n + 1));
    std::vector<double> X(n, 0);
    
    //std::cout << "Введите матрицу A: " << std::endl;

    for( int i = 0; i < n; i++)
        for( int j = 0; j < n; j++)
            std::cin >> A[i][j];

    //std::cout << "Введите вектор правых частей: ";
    
    for(int i = 0; i < n; i++ )
        std::cin >> A[i][n];

    //std::cout << "Введите точность: ";
    std::cin >> eps;

    std::cout << "Число итераций метода простых итераций: " << IterMethod(A, X, eps) << std::endl;
    for(int i = 0; i < n; i++)
        std::cout << "  x_" << i + 1 << " = " << X[i] << std::endl;    

    return 0;
}
