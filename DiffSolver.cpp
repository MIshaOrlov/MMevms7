#include "DiffSolver.h"

// Правая часть уравнения: e^(x^2)
double ex2(double x) {
    return exp(x * x);
}


// Основная функция решения методом конечных разностей
int SolveFiniteDifference(int N, double* A, double* B, double* u, double* x, double h, double* L, double* U){
    // Инициализация сетки x
    for (int i = 0; i < N; ++i) {
        x[i] = i * h;
        //std::cout << x[i] << " ";
    }
    //std::cout << std::endl;

    // Инициализация решения и правой части
    for (int i = 0; i < N; ++i) {
        u[i] = 0.0;            // Начальные условия для u
        B[i] = ex2(x[i]);      // Правая часть системы
    }

    // Инициализация матрицы A (нулевая матрица)
    for (int i = 0; i < N * N; ++i) {
        A[i] = 0.0;
    }

    // Заполняем матрицу A согласно разностной схеме для производных
    for (int i = 2; i < N - 2; ++i) {
        A[i * N + (i - 2)] = 1.0 / (h * h * h * h) - 1.0 / (2 * h * h * h);
        A[i * N + (i - 1)] = -4.0 / (h * h * h * h) + 2.0 / (2 * h * h * h);
        A[i * N + i] = 6.0 / (h * h * h * h) + sin(x[i]);
        A[i * N + (i + 1)] = -4.0 / (h * h * h * h) - 2.0 / (2 * h * h * h);
        A[i * N + (i + 2)] = 1.0 / (h * h * h * h) + 1.0 / (2 * h * h * h);
    }

    // Граничные условия
    
    // u(0) = 0, u'(0) = 0
    A[0 * N + 0] = 1.0;
    A[1 * N + 0] = 1.0 / h;
    A[1 * N + 1] = -1.0 / h;

    B[0] = 0.0;
    B[1] = 0.0;

    // u(1) = 1, u'(1) = 0
    A[(N - 1) * N + (N - 1)] = 1.0;
    A[(N - 2) * N + (N - 1)] = 1.0 / h;
    A[(N - 2) * N + (N - 2)] = -1.0 / h;

    B[N - 1] = 1.0;
    B[N - 2] = 0.0;



    LU_decomposition(N, A,L,U);
    solve_LU(N, L, U, u, B);

    

    return 1;
}
