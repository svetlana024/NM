#include <vector>
#include <math.h>
#include <iostream>

void PrintF1(std::vector<double>& F1)
{
    std::cout << "Приближающий многочлен 1-й степени F1(x) = ";
    if(F1[0] != 0)
        std::cout << F1[0];
      
    if(F1[1] < 0)
        std::cout << " - " << -F1[1] << "x" << std::endl;
    else if(F1[1] > 0)
        std::cout << " + " << F1[1] << "x" << std::endl;
    else std::cout << std::endl;
}

void PrintF2(std::vector<double>& F2)
{
    std::cout << "Приближающий многочлен 2-й степени F2(x) = ";
    if(F2[0] != 0)
        std::cout << F2[0];
      
    if(F2[1] < 0)
        std::cout << " - " << -F2[1] << "x";
    else if(F2[1] > 0)
        std::cout << " + " << F2[1] << "x";

    if(F2[2] < 0)
        std::cout << " - " << -F2[2] << "x^2" << std::endl;
    else if(F2[2] > 0)
        std::cout << " + " << F2[2] << "x^2" << std::endl;
    else std::cout << std::endl;
}

void Swap(std::vector<std::vector<double>>& Matrix, int i, int k)
{
    std::vector<double> H = Matrix[i];
    Matrix[i] = Matrix[k];
    Matrix[k] = H;
}

void Back(std::vector<std::vector<double>>& Matrix, std::vector<double>& x, int n)
{
    double sum = 0;
    int j = 0;

    for( int i = n -1; i >= 0; i--)
    {
        j = i + 1;
        while(j < n){
            sum += x[j] * Matrix[i][j];
            j++;
        }
        x[i] = ( Matrix[i][n] - sum) / Matrix[i][i];
        sum = 0;
    }   
}
        
void Straight(std::vector<std::vector<double>>& Matrix, int n)
{
    int k = 0;
    double Kaf = 0;

    for( int i = 0; i < n - 1; i++)
    {
        k = i;
        while( Matrix[k][i] == 0 )
            k++;
        Swap(Matrix, i, k);
        for( int j = i + 1; j < n; j++)
        {
            Kaf = Matrix[j][i] / Matrix[i][i];
            for( int h = i; h <= n; h++)
                Matrix[j][h] -= Kaf * Matrix[i][h];
        }
    }
}

int main()
{
    int n;
    double mistake = 0;

    //std::cout << "Введите количество интерполяционных узлов: ";

    std::cin >> n;
    std::vector<double> x(n);
    std::vector<double> y(n);
    std::vector<double> F1(2, 0);
    std::vector<double> F2(3, 0);
    std::vector<std::vector<double>> A1(2, std::vector<double> (3, 0)); 
    std::vector<std::vector<double>> A2(3, std::vector<double> (4, 0));

    A1[0][0] = n;
    A2[0][0] = n;

   //std::cout << "Введите множество точек - интерполяционных узлов:";
    for(int i = 0; i < n; i++)
    {
        std::cin >> x[i];
        A1[0][1] += x[i];
        A1[1][1] += pow(x[i],2);
        A2[1][2] += pow(x[i], 3);
        A2[2][2] += pow(x[i], 4);
    }
    A1[1][0] = A1[0][1];
    A2[0][1] = A1[0][1];
    A2[1][0] = A2[0][1];
    A2[0][2] = A1[1][1];
    A2[1][1] = A2[0][2];
    A2[2][0] = A2[0][2];
    A2[2][1] = A2[1][2];
    
    //std::cout << "Введите множество значений функции в узлах интерполяции:";
    for(int i = 0; i < n; i++)
    {
        std::cin >> y[i];
        A1[0][2] += y[i];
    }
    
    A2[0][3] = A1[0][2];
    for(int i = 0; i < n; i++)
    {
        A1[1][2] += x[i]*y[i];
        A2[2][3] += y[i]*pow(x[i], 2);
    }
    A2[1][3] = A1[1][2];

    Straight(A1, 2);
    Back(A1, F1, 2);
    PrintF1(F1);

    for(int i = 0; i < n; i++)
        std::cout << "F1(x_" << i + 1 << ") = " << F1[0] + F1[1]*x[i] << std::endl;

    for(int i = 0; i < n; i++)
        mistake += pow(F1[0] + F1[1]*x[i] - y[i], 2);
    std::cout << "Сумма квадратов ошибок Ф = " << mistake << std::endl;
    mistake = 0;

    Straight(A2, 3);
    Back(A2, F2, 3);
    PrintF2(F2);

    for(int i = 0; i < n; i++)
        std::cout << "F2(x_" << i + 1 << ") = " << F2[0] + F2[1]*x[i] + F2[2]*pow(x[i], 2) << std::endl;

    for(int i = 0; i < n; i++)
        mistake += pow(F2[0] + F2[1]*x[i] + F2[2]*pow(x[i], 2) - y[i], 2);
    std::cout << "Сумма квадратов ошибок Ф = " << mistake << std::endl;


    return 0;
}

