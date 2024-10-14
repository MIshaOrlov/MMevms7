#ifndef MatrixHandlers_hpp
#define MatrixHandlers_hpp


#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <time.h>

// Функция для вычисления элемента матрицы по формуле
double calculateElement(int k, int n, int i, int j);
// Функция для создания по формуле
void CreateMatrix(double* matrix, int n, int k);


// Функция для чтения из файла
int inputMatrix(double* matrix, int n, const std::string& filename);


// Функция для освобождения памяти, выделенной для матрицы
void deleteMatrix(double* matrix);

// Функция для вывода матрицы
void printMatrix(double* matrix, int n,int m);

// Функция для вывода вектора
void printVector(double* matrix, int n);

// Функция для транспонирования матрицы
void transposeMatrix(double* matrix, int n);

// Функция для вычисления максимальной столбичной нормы матрицы
double calculateColumnNorm(double* matrix, int n);

// Функция для вычисления максимальной строчной нормы матрицы
double calculateRowNorm(double* matrix, int n);
#endif /* MatrixHandlers_hpp */
