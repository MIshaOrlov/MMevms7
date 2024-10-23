#include "DiffSolver.h"

// Правая часть уравнения: e^(x^2)
double ex2(double x) {
    return exp(x * x);
}

// Основная функция решения методом конечных разностей
int SolveFiniteDifference(int N, double* A, double* B, double* u, double* x, double h, double* adjoint) {
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

    // Вывод матрицы A
    //printMatrix(A, N, N);

    // Инверсия матрицы методом отражений
    int res = ReflectionInverse(N, A, adjoint);

    if (res == -1) {
        std::cout << "Не существует обратной матрицы" << std::endl;
        return -1;
    }

    // Умножение инвертированной матрицы на вектор B для получения u
    multiplyMatrixByVector(N, adjoint, B, u);

    // Вывод результатов в файл
   /* std::ofstream outfile("solution.txt");
    for (int i = 0; i < N; ++i) {
        outfile << x[i] << " " << u[i] << std::endl;
    }
    outfile.close();

    std::cout << "Решение сохранено в файл 'solution.txt'." << std::endl;*/

    return 1;
}
