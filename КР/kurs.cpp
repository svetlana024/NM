#include <iostream>
#include <math.h>
#include <vector>
#include <iterator>
#include "matrix.h"

void Housholder( Matrix* R, Matrix* Q)
{
    int n = R-> n;
    int m = R-> m;
    ClearMatrix(Q);
    InitE(Q, n);
    for( int i = 0; i < m - 1; i++)
    {
        std::vector<double> v(n, 0);
        double sum = 0;
        for(int j = R-> points[i]; j < R-> points[i + 1]; j++)
        {
            if( R-> columns[j] == i)
            {
                sum += pow( R-> values[j], 2);
                v[i] = R-> values[j];
            }
            else if( R-> columns[j] == i + 1)
            {
                sum += pow( R-> values[j], 2);
                v[i] += Sign(R-> values[j-1])*pow(sum, 0.5);
                v[i + 1] = R-> values[j];
            }
        }
        Matrix Housholder;

        InitHousholder(v, &Housholder);
        TransMatrix(R);

        MultMatrix(&Housholder, R);

        MultMatrix(Q, &Housholder);
        InitMatrixWithMatrix(Q, &Housholder);

        if( i == m - 2)
        {
            Matrix A_new;
            MultMatrix(R, Q, &A_new);
            InitMatrixWithMatrix(R, &A_new);
        }
        TransMatrix(R);
    }
}

Matrix QR_algorithm( Matrix* A, std::vector<Complex>& Lyambda)
{
    Matrix Q;

    int it_count = 0;
    double kr = 1;
    double eps = 1e-6;
    int n = A-> m;

    //std::vector<Complex> Lyambda(n, {0, 0});
    TransMatrix(A);

    while(kr > eps)
    {
        Housholder(A, &Q);
        for( int i = 0; i < n; i++ )
        {
            kr = 0;
            for(int j = A-> points[i]; j < A-> points[i + 1]; j++)
                if( A-> columns[j] > i)
                    kr += pow(A-> values[j], 2);
            kr = pow(kr, 0.5);

            if(kr <= eps)
            {
                Lyambda[i].real = 0;
                Lyambda[i].imag = 0;
                
                for( int k = A-> points[i]; k < A-> points[i + 1]; k++)
                    if( A-> columns[k] == i )
                    {
                        Lyambda[i].real = A-> values[k];
                        break;
                    }
                //if( i == 1)
                //    it_count += 500;
            }
            else break;
            /*else if( i == 0)
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
            }*/
        }
        it_count ++;
    }
    //if( it_count > 501)
    //    it_count -= 500;
    TransMatrix(A);
    std::cout << "Верхняя треугольная матрица A в подпространстве Крылова: " << std::endl;
    PrintMatrix(A);
    std::cout << std::endl;

    //PrintMatrix(&Q);
    //std::cout << "Число итераций: " << it_count << std::endl;
    std::cout << "Собственные значения матрицы A:" << std::endl;
    for(int i = 0; i < n; i++)
        if(Lyambda[i].imag == 0)
            std::cout << std::setw(10) << std::right << "lambda_" << i << " = " << Lyambda[i].real << std::endl;
        else
        {
            std::cout << std::setw(10) << std::right << "lambda_" << i << " = " << Lyambda[i].real;
            std::cout<< std::showpos << Lyambda[i].imag << "i" << std::endl;
        }
                   
    return Q;

}

void Check( Matrix* A, std::vector<Complex>& Lyambda)
{
    int n = Lyambda.size();
    double trace = 0;

    std::cout << "След матрицы H: ";
    for( int i = 0; i < n; i++)
    {
        for( int j = A-> points[i]; j < A-> points[i + 1]; j++)
            if( A-> columns[j] == i)
            {
                std::cout << std::showpos << A-> values[j] << " ";
                trace += A-> values[j];
            }
    }

    std::cout << " = " << trace << std::endl;
    trace = 0;
    std::cout << "Сумма собственных значений : ";
    for( int i = 0; i < n; i++)
    {
        std::cout << std::showpos << Lyambda[i].real << Lyambda[i].imag << " ";
        trace += Lyambda[i].real;
    }

    std::cout << " = " << trace << std::endl;
}


void ArnoldiIteration(Matrix* A)
{
    int dim = 3, points = 0;
    double eps = 1e-10;

    Matrix Q, H;
    std::vector<double> q(A-> m, 0);
    q[0] = 0.78289;
    q[1] = 0.62216;

    SetDim(&Q, dim + 1, A-> n);
    SetDim(&H, dim, dim + 1);
    SetRow(&Q, q);
    
    for( int i = 0; i < dim; i++)
    {
        std::vector<double> v;
        MultMatVec(A, q, v);

        for( int j = 0; j < i + 1; j++)
        {
            double element = 0;
            for( int k = Q.points[j]; k < Q.points[j + 1]; k++)
                element += Q.values[k] * v[Q.columns[k]];
            if(element != 0){
                H.values.push_back(element);
                H.columns.push_back(j);
                points++;
            }
            for( int k = Q.points[j]; k < Q.points[j + 1]; k++)
                v[Q.columns[k]] -= element * Q.values[k];
        }

        double norma = Norma(v);
        if( norma != 0)
        {
            H.values.push_back(norma);
            H.columns.push_back(i+1);
            points++;
        }
        H.points.push_back(points);
        if( norma > eps)
        {
            Normalization(v, norma, q);
            SetRow(&Q, q);
        }
        else break;
    }
    Matrix TQ, AA;
    InitMatrixWithMatrix(&TQ, &Q);
    RemRow(&TQ);
    TransMatrix(&TQ);

    TransMatrix(&Q);        //Basis m(n+1)
    TransMatrix(&H, &AA);
    TransMatrix(&H);        //Hessenberg (n+1)n

    std::cout << "Матрица Хессенберга: " << std::endl;
    PrintMatrix(&H);
    std::cout << "Базис Крылова: " << std::endl;
    PrintMatrix(&Q);

    RemRow(&H);
    RemRow(&AA);
    std::vector<Complex> Lyambda(H.m, {0, 0});
    Matrix SV = QR_algorithm(&H, Lyambda);      //Nortmal Hessenberg matrix with removed last row

    std::cout << "Собственные векторы: " << std::endl;
    
    /*for( int i = SV.m; i < A-> m; i++)
        SV.points.push_back(SV.points[SV.m]);
    SV.m = A-> m;
    MultMatrix(&Q, &SV);
    PrintMatrix(&SV);*/

    //RemRow(&TQ);
    //PrintMatrix(&TQ);
    MultMatrix(&TQ, &SV);
    //std::cout << "Собственные векторы после преобразования подобия: " << std::endl;
   
    PrintMatrix(&SV);
    std::cout << "Проверка: " << std::endl;
    Check(&AA, Lyambda);
}


int main(void)
{

    Matrix matrix;
    InitMatrix(&matrix);
    ArnoldiIteration(&matrix);
    ClearMatrix(&matrix);

    return 0;
}

    
