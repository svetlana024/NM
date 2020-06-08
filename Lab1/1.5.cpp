#include <iostream>
#include <vector>
#include <cmath>

typedef struct Complex{
    double real;
    double imag;
} Complex;

void PrintMatrix(std::vector<std::vector<double>>& Matrix)
{
    for(int i = 0; i < Matrix.size(); i++){
        for(int j = 0; j < Matrix[0].size(); j++)
            std::cout << Matrix[i][j] <<  " ";
        std::cout << std::endl;
    }
}

void MultMatrix(std::vector<std::vector<double>>& F, std::vector<std::vector<double>>& S)
{
    int n = F.size();
    int m = S[0].size();
    int l = S.size();
    double a = 0;
    std::vector<double> H(n);

    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
        {
            for(int k = 0; k < l; k++)
                a += F[j][k]*S[k][i];
            H[j] = a;
            a = 0;
        }
        for(int j = 0; j < n; j++)
            S[j][i] = H[j];
    }
}

void MultMatrix(std::vector<std::vector<double>>& F, std::vector<std::vector<double>>& S, std::vector<std::vector<double>>& R)
{
    int n = F.size();
    int m = S.size();
    int l = S[0].size();
    int k = R.size();
    int h = R[0].size();
    double a = 0;

    if(n - k > 0)
        for(int i = 0; i < n - k; i++)
            R.insert(R.end(), std::vector<double>(l, 0));
    if(l - h > 0)
        for(int i = 0; i < n; i++)
            for(int j = 0; j < l - h; j++)
                R[i].insert(R[i].end(), 0);

    for(int i = 0; i < l; i++)
        for(int j = 0; j < n; j++)
        {
            for(int k = 0; k < m; k++)
                a += F[j][k]*S[k][i];
            R[j][i] = a;          
            a = 0;
        }
    if(k - n > 0)
        R.erase(R.end()- k + n, R.end());
    if(h - l > 0)
        for(int i = 0; i < n; i++)
            R[i].erase(R[i].end() - h + l, R[i].end());
}

int Sign(double a)
{
    return
        (a < 0) ? -1
        : (a > 0) ? 1
        : 0;
}
  
void Transporant(std::vector<std::vector<double>>& M)
{
    double h = 0;
    int n = M.size();

    for(int i = 0; i < n; i++)
        for(int j = i + 1; j < n; j++)
        {
            h = M[i][j];
            M[i][j] = M[j][i];
            M[j][i] = h;
        } 
} 

void HouseHolder(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& Q)
{
    int n = A.size();
    double sum = 0;
    double kaf = 0;
    std::vector<std::vector<double>> Nju(n, std::vector<double>(1, 0));
    std::vector<std::vector<double>> NjuTransp(1, std::vector<double>(n, 0));
    std::vector<std::vector<double>> H(n, std::vector<double>(n, 0));
    std::vector<std::vector<double>> R(1, std::vector<double>(1, 0));
    
    for(int i = 0; i < n - 1; i ++)
    {
        for(int j = 0; j < i; j++)
        {
            Nju[j][0] = 0;
            NjuTransp[0][j] = 0;
        }
        for(int j = i; j < n; j++)
        {
            Nju[j][0] = A[j][i];
            NjuTransp[0][j] = Nju[j][0];
            sum += pow(A[j][i], 2);
        }

        Nju[i][0] += Sign(Nju[i][0])*pow(sum, 0.5);
        NjuTransp[0][i] = Nju[i][0];
        MultMatrix(NjuTransp, Nju, R);
        kaf = -2/R[0][0];
        MultMatrix(Nju, NjuTransp, R);
    
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                if(i == j)
                    H[i][j] = 1 + R[i][j]*kaf;
                else
                    H[i][j] = R[i][j]*kaf;
        MultMatrix(H, A);
        if(i == 0)
            Q = H;
        else
            MultMatrix(Q, H);
        Q = H;
        sum = 0;
    }
    H = A;
    MultMatrix(H, Q, A);     
}

void GetComplexSV(std::vector<std::vector<double>>& A, std::vector<Complex>& Lyambda, int i)
{
    double Disc = pow(A[i+1][i+1] + A[i][i], 2) - 4*(A[i][i]*A[i+1][i+1] - A[i][i+1]*A[i+1][i]);

    Lyambda[i].real = (A[i+1][i+1] + A[i][i])/2;
    Lyambda[i].imag = pow(fabs(Disc), 0.5)/2;
    Lyambda[i+1].real = Lyambda[i].real;
    Lyambda[i+1].imag = -Lyambda[i].imag;
}

int main()
{
    int n, it_count = 0;
    double eps, kr = 1;

    //std::cout << "Введите размерность матрицы: ";
    std::cin >> n;

    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    std::vector<std::vector<double>> Q(n, std::vector<double>(n));
    std::vector<Complex> Lyambda(n, {0, 0});
    
    //std::cout << "Введите матрицу A: " << std::endl;

    for( int i = 0; i < n; i++)
        for( int j = 0; j < n; j++)
            std::cin >> A[i][j];

    //std::cout << "Введите точность: " << std::endl;
    std::cin >> eps;  

    while(kr > eps)
    {
        HouseHolder(A, Q);
        for(int i = 0; i < n; i++ )
        {
            kr = 0;
            for(int j = i + 1; j < n; j++)
                kr += pow(A[j][i], 2);
            kr = pow(kr, 0.5);
            if(kr <= eps)
                Lyambda[i].real = A[i][i];
            else if( i == 0)
                break;
            else if(i != n - 1)
            {
                double real = Lyambda[i].real;
                double img = Lyambda[i].imag;
                GetComplexSV(A, Lyambda, i);
                kr = pow(pow(real - Lyambda[i].real, 2) + pow(img - Lyambda[i].imag, 2), 0.5);
                i++;
                if(kr > eps)
                    break;
            }
        }
        it_count ++;
    }
    //std::cout << "Число итераций: " << it_count << std::endl;
    std::cout << "Собственные значения матрицы A:" << std::endl;
    for(int i = 0; i < n; i++)
        if(Lyambda[i].imag == 0)
            std::cout << "lyambda_" << i << " = " << Lyambda[i].real << std::endl;
        else
            std::cout << "lyambda_" << i << " = " << Lyambda[i].real << " + " << Lyambda[i].imag << "i" << std::endl;
                   
    return 0;
}
