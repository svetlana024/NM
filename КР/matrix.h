#ifndef MATRIX_H
#define MATRIX_H

#include <iomanip>

typedef struct
{
    int m;
    int n;
    std::vector<double> values;
    std::vector<int> columns;
    std::vector<int> points;
} Matrix;

typedef struct Complex{
    double real;
    double imag;
} Complex;

void InitHousholder(std::vector<double>& v, Matrix* matrix){

    int n = v.size();
    matrix-> m = n;
    matrix-> n = n;

    int points = 0;
    double kaf = 0;
    double element;
    matrix-> points.push_back(points);

    for( int i = 0; i < n; i++)
        kaf += v[i] * v[i];

    for( int i = 0; i < n; i++)
    {
        for( int j = 0; j < n; j++)
        {
            if( i == j)
                element = 1 - 2 * v[i]*v[j]/kaf;
            else
                element = -2 * v[i]*v[j]/kaf;

            if( element != 0)
            {
                matrix-> values.push_back(element);
                matrix-> columns.push_back(j);
                points++;
            }
        }
        matrix-> points.push_back(points);
    }
}

void InitE( Matrix* E, int n)
{
    E-> m = n;
    E-> n = n;

    for( int i = 0; i < n; i++)
    {   
        E-> values.push_back(1);
        E-> columns.push_back(i);
        E-> points.push_back(i);
    }
    E-> points.push_back(n);
}

int Sign( double value)
{
    return (value > 0) ? 1 : -1;
}

void ClearMatrix( Matrix* matrix)
{
    matrix-> values.clear();
    matrix-> columns.clear();
    matrix-> points.clear();
}

void InitMatrix( Matrix* matrix)
{
    int value;
    int count = 0;
    int m;

    std::cout << "Введите размерность матрицы: ";
    std::cin >> m;

    matrix-> points.push_back(0);
    matrix-> m = m;
    matrix-> n = m;

    std::cout << "Введите матрицу: " << std::endl;
    for( int i = 0; i < m; i++)
    {
        for( int j = 0; j < m; j++)
        {
            std::cin >> value;
            if( value != 0)
            {
                matrix-> values.push_back(value);
                matrix-> columns.push_back(j);
                count++;
            }
        }
        matrix-> points.push_back(count);
    }
}

void InitMatrixWithMatrix( Matrix* m1, Matrix* m2)
{
    m1-> m = m2-> m;
    m1-> n = m2-> n;
    m1-> values = m2-> values;
    m1-> columns = m2-> columns;
    m1-> points = m2-> points;
}


void TransMatrix( Matrix* m)
{
    int points = 0;
    int n = m-> n;

    Matrix trans;
    trans.points.push_back(points);
    trans.m = m-> n;
    trans.n = m-> m;

    for( int i = 0; i < n; i++)
    {
        int str = 0;
        for( int j = 0; j < m-> points.size() - 1; j++)
        {
            for( int k = m-> points[j]; k < m-> points[j + 1]; k++)
            {
                if( m-> columns[k] == i)
                {
                    trans.values.push_back(m-> values[k]);
                    trans.columns.push_back(str);
                    points++;
                }
            }
            str++;
        }
        trans.points.push_back(points);
    }
    ClearMatrix(m);
    m-> m = trans.m;
    m-> n = trans.n;
    m-> values = trans.values;
    m-> columns = trans.columns;
    m-> points = trans.points;
}


void TransMatrix( Matrix* m, Matrix* trans)
{
    int points = 0;
    int n = m-> n;

    trans-> points.push_back(points);
    trans-> m = m-> n;
    trans-> n = m-> m;

    for( int i = 0; i < n; i++)
    {
        int str = 0;
        for( int j = 0; j < m-> points.size() - 1; j++)
        {
            for( int k = m-> points[j]; k < m-> points[j + 1]; k++)
            {
                if( m-> columns[k] == i)
                {
                    trans-> values.push_back(m-> values[k]);
                    trans-> columns.push_back(str);
                    points++;
                }
            }
            str++;
        }
        trans-> points.push_back(points);
    }
}

void PrintMatrix( Matrix* matrix)
{
    int id = 0;
    for( int i = 0; i < matrix-> points.size()-1; i++)
    {
        for( int j = matrix-> points[i]; j < matrix-> points[i+1]; j++)
        {
            for( int k = id; k < matrix-> columns[j]; k++)
                std::cout << std::setw(15) << std::right << "0";
            std::cout << std::setw(15) << std::right << matrix-> values[j];
            id = matrix-> columns[j] + 1;
        }

        for( int j = id; j < matrix-> n; j++)
            std::cout << std::setw(15) << std::right << "0";
        id = 0;
        std::cout << std::endl;
    }
}

void MultMatVec( Matrix* m, std::vector<double>& v, std::vector<double>& res)
{
    double element = 0;
    int points = 0;
    
    for( int i = 0; i < m-> points.size() - 1; i++)
    {
        for( int k = m-> points[i]; k < m-> points[i + 1]; k++)
                element += m-> values[k] * v[m-> columns[k]];
        res.push_back(element);
        element = 0;
    }
}

void MultMatrix( Matrix* m1, Matrix* m2, Matrix* res)
{
    double element = 0;
    int points = 0;

    Matrix trans;
    TransMatrix(m2, &trans);

    res-> points.push_back(points);
    res-> m = m1-> m;
    res-> n = trans.m; 

    for( int i = 0; i < m1-> points.size() - 1; i++)
    {
        for( int j = 0; j < trans.points.size() - 1; j++)
        {
            for( int k = m1-> points[i]; k < m1-> points[i + 1]; k++)
            {
                for( int l = trans.points[j]; l < trans.points[j + 1]; l++)
                {
                    if( m1-> columns[k] == trans.columns[l])
                        element += m1-> values[k] * trans.values[l];
                }
            }
            if( element != 0)
            {
                res-> values.push_back(element);
                res-> columns.push_back(j);
                points++;
                element = 0;
            }
        }
        res-> points.push_back(points);
    }
}


void MultMatrix( Matrix* m1, Matrix* m2)
{
    double element = 0;
    int points = 0;

    Matrix trans;
    TransMatrix(m2, &trans);

    ClearMatrix(m2);
    m2-> m = m1-> m;
    m2-> points.push_back(points);

    for( int i = 0; i < m1-> points.size() - 1; i++)
    {
        for( int j = 0; j < trans.points.size() - 1; j++)
        {
            for( int k = m1-> points[i]; k < m1-> points[i + 1]; k++)
            {
                for( int l = trans.points[j]; l < trans.points[j + 1]; l++)
                {
                    if( m1-> columns[k] == trans.columns[l])
                        element += m1-> values[k] * trans.values[l];
                }
            }
            if( element != 0)
            {
                m2-> values.push_back(element);
                m2-> columns.push_back(j);
                points++;
                element = 0;
            }
        }
        m2-> points.push_back(points);
    }
}

void SetDim(Matrix* matrix, int m, int n)
{
    matrix-> m = m;
    matrix-> n = n;
    matrix-> points.push_back(0);
}

void RemRow( Matrix* matrix)
{
    matrix-> m = matrix-> m - 1;
    int shift = matrix-> points[matrix-> points.size() - 1] - matrix-> points[matrix-> points.size() - 2];
    matrix-> points.pop_back();
    
    while( shift > 0)
    {
        matrix-> values.pop_back();
        matrix-> columns.pop_back();
        shift--;
    }
}

void SetRow( Matrix* matrix, std::vector<double>& v)
{
    int n = v.size();
    int points = matrix-> points[matrix-> points.size() - 1];
    for( int i = 0; i < n; i++)
    {
        if(v[i] != 0)
        {
            matrix-> values.push_back(v[i]);
            matrix-> columns.push_back(i);
            points++;
        }
    }
    matrix-> points.push_back(points);
}

double Norma( std::vector<double>& v)
{
    double norma = 0;
    for( int i = 0; i < v.size(); i++)
        norma += pow(v[i], 2);
    return pow(norma, 0.5);
}

void Normalization( std::vector<double>& from, double norma, std::vector<double>& to)
{
    for(int i = 0; i < from.size(); i++)
        to[i] = from[i] / norma;
}

void Print( Matrix* matrix)
{
    std::vector<int>::iterator it;
    std::vector<double>::iterator it_v;
    std::cout << "m x n: " << matrix-> m << " " << matrix-> n << std::endl;
    std::cout << "values: ";
    for(it_v = matrix-> values.begin(); it_v != matrix-> values.end(); it_v++)
        std::cout << *it_v << " ";
    std::cout << std::endl << "columns: ";
    for(it = matrix-> columns.begin(); it != matrix-> columns.end(); it++)
        std::cout << *it << " ";
    std::cout << std::endl << "points: ";
    for(it = matrix-> points.begin(); it != matrix-> points.end(); it++)
        std::cout << *it << " ";
}


void GetComplexSV( Matrix* A, std::vector<Complex>& Lyambda, int i)
{
    double Disc = pow(A-> values[A-> points[i + 1] + i + 1] + A-> values[A-> points[i] + i], 2) - 4*(A-> values[A-> points[i] + i] * A-> values[A-> points[i + 1] + i + 1] - A-> values[A-> points[i] + i + 1] * A-> values[A-> points[i + 1] + i]);

    Lyambda[i].real = (A-> values[A-> points[i + 1] + i + 1] + A-> values[A-> points[i] + i])/2;
    Lyambda[i].imag = pow(fabs(Disc), 0.5)/2;
    Lyambda[i+1].real = Lyambda[i].real;
    Lyambda[i+1].imag = -Lyambda[i].imag;
}

#endif
