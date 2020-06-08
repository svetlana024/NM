#include <iostream>
#include <vector>

void PrintMatrix(std::vector<std::vector<double>>& Matrix)
{
    int n = Matrix.size();
    int m = Matrix[0].size();
    for( int i = 0; i < n; i++)
    {
        for( int j = 0; j < m; j++)
            std::cout << Matrix[i][j] <<  " ";
        std::cout << std::endl;
    }
}

void SwapStrings(std::vector<std::vector<double>>& Matrix, int i, int k)
{
    int n = Matrix.size();
    std::vector<double> H = Matrix[i];
    Matrix[i] = Matrix[k];
    Matrix[i] = H;
}

void ChangeWithMax(std::vector<std::vector<double>>& Matrix, int i)
{
    int n = Matrix.size();
    double max = Matrix[i][i]*Matrix[i][i];
    int max_id = i;

    for(int j = i + 1; j < n; j++)
    {
        if( Matrix[j][i]*Matrix[j][i] > max)
        {
            max = Matrix[j][i]*Matrix[j][i];
            max_id = j;
        }
    }
    if(max_id != i)
        SwapStrings(Matrix, i, max_id);
}

void FindX(std::vector<std::vector<double>>& Matrix, std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& b)
{
    int n = Matrix.size();
    int m = x[0].size();
    int j;
    
    for(int k = 0; k < m; k++)
    {
        for(int i = 0; i < n; i++)                                  //решение уравнения Ly = b ---> y записываем в b
        {
            j = i - 1;
            while(j >= 0)
            {
                b[i][k] -= Matrix[i][j] * b[j][k];
                j--;  
            }
        }
    }
    for(int k = 0; k < m; k++)
    {
        for( int i = n - 1; i >= 0; i--)                                //решение уравнения Ux = b обратным ходом
        {
            j = i + 1;
            while(j < n){
                b[i][k] -= Matrix[i][j] * x[j][k];
                j++;
            }
            x[i][k] = b[i][k] / Matrix[i][i];
        }
    }
}
                
void LU(std::vector<std::vector<double>>& Matrix)
{
    int n = Matrix.size();
    int main_id = 0;
    PrintMatrix(Matrix);
    ChangeWithMax(Matrix, 0);
    
    for( int j = 1; j < n; j++)                                     // i = 0
        Matrix[j][0] = Matrix[j][0]/Matrix[0][0];
    
    for( int i = 1; i < n; i++)                                     // i = 1..n
    {
        for( int j = i; j < n; j++)                                 // u[i][j]
        {
            for( int k = 0; k < i; k++)
                Matrix[i][j] -= Matrix[i][k]*Matrix[k][j];
        }
        ChangeWithMax(Matrix, i);

        for( int j = i + 1; j < n; j++)                                 //l[j][i]
        {
            for( int k = 0; k < i; k++)
                Matrix[j][i] -= Matrix[j][k] * Matrix[k][i];
            Matrix[j][i] /= Matrix[i][i];
        }
    }
}
    
double Determinant(std::vector<std::vector<double>>& Matrix)
{
    int det = 1;
    for( int i = 0; i < Matrix.size(); i++)
        det *= Matrix[i][i];
    return det;
}

void L(std::vector<std::vector<double>>& Matrix)
{
    int n = Matrix.size();  
    int count = 0;

    for( int i = 0; i < n; i++)
    {
        count = i;
        for( int j = 0; j < n; j++)
        {
            if( count > 0)
            {
                std::cout << Matrix[i][j] << " ";
                count--;
            }
            else if( count == 0)
            {
                std::cout << "1 ";
                count--;
            }
            else
               std::cout << "0 ";
        }
        std::cout << std::endl;
    }
}

void U(std::vector<std::vector<double>>& Matrix)
{
    int n = Matrix.size();
    int count = 0;

    for( int i = 0; i < n; i++)
    {
        count = i;
        for( int j = 0; j < n; j++)
        {
            if( count > 0)
            {
                std::cout << "0 ";
                count--;
            }
            else
                std::cout << Matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main()
{
    int n;
    double element;

    //std::cout << "Введите размерность матрицы: ";
    std::cin >> n;

    std::vector<std::vector<double>> A(n, std::vector<double>(n + 1));
    std::vector<std::vector<double>> x(n, std::vector<double>(1));
    std::vector<std::vector<double>> b(n, std::vector<double>(1));
    std::vector<std::vector<double>> ReversedA(n, std::vector<double>(n, 0));
    std::vector<std::vector<double>> ZeroMatrix(n, std::vector<double>(n ,0));
    
    //std::cout << "Введите матрицу A: " << std::endl;

    for( int i = 0; i < n; i++)
        for( int j = 0; j < n; j++)
        {
            std::cin >> A[i][j];
            if( i == j)
                ZeroMatrix[i][j] = 1;
        }

    //std::cout << "Введите результирующий столбец: ";
    
    for(int i = 0; i < n; i++ )
        std::cin >> A[i][n];
    LU(A);
    std::cout << "LU-разложение:" << std::endl;
    std::cout << "Нижняя треугольная матрица L:" << std::endl;
    L(A);
    std::cout << std::endl;
    std::cout << "Верхняя треугольная матрица U:" << std::endl;
    U(A);
    std::cout << std::endl;
    std::cout << "Решение системы с помощью LU-разложения: " << std::endl;
    for( int i = 0; i < n; i++)
        b[i][0] = A[i][n];
    FindX(A, x, b);
    
    for(int i = 0; i < n; i++)
        std::cout << "x_" << i + 1 << " = " << x[i][0] << std::endl;

    std::cout << "Детерминант: " << Determinant(A) << std::endl;
    std::cout << std::endl;

    FindX(A, ReversedA, ZeroMatrix);
    std::cout << "Обратная матрица: " << std::endl;
    PrintMatrix(ReversedA);
    
    return 0;
}
