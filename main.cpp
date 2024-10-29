/*
 3 Построить разностную схему со вторым порядком аппроксимации и найти ее решение при раз-
 личных значениях h:
 u(4) + u(3)+ sinx · u = ex2,
 u(0) = u′(0) = 0, u(1) = 1, u′(1) = 0
 Исследовать построенную разностную схему на устойчивость и сходимость.
 */


#include <iostream>
#include <cstdlib>        // Для std::stod
#include <ctime>          // Для измерения времени работы алгоритма
#include <fstream>        // Для сохранения результатов в файл
#include "Matrix.h"       // Библиотека для работы с матрицами
#include "ReflectionMethod.h" // Для метода отражений
#include "DiffSolver.h"   // Для метода конечных разностей

using namespace std;

int main(int argc, char* argv[]) {
    // Проверяем, что передано два аргумента: шаг h и имя выходного файла
    if (argc < 3) {
        cout << "Недостаточно аргументов командной строки." << endl;
        cout << "Использование: " << argv[0] << " <h> <output_file>" << endl;
        return 1;
    }

    // Преобразуем переданный аргумент h в число с плавающей запятой
    double h = stod(argv[1]);

    // Имя файла для сохранения решения
    const char* output_file = argv[2];

    // Определяем количество узлов на сетке
    int n = static_cast<int>(1.0 / h) + 1;

    // Вывод аргументов
    cout << "Полученные аргументы: h = " << h << ", n = " << n << endl;

    // Выделяем память под матрицы и векторы
    double* matrix = new double[n * n];   // Матрица A
    double* b = new double[n];            // Вектор правой части
    double* x = new double[n];            // Вектор сетки (x)
    double* u = new double[n];            // Вектор решения u
    // Массивы для коэффициентов метода прогонки
    double* L = new double[n * n];
    double* U = new double[n * n];
    
    
    
    // Измеряем время выполнения
    clock_t t;
    t = clock();

    // Решаем задачу методом конечных разностей
    int res = SolveFiniteDifference(n, matrix, b, u, x, h, L,U);

    // Если не удалось найти решение (обратную матрицу), выводим сообщение
    if (res == -1) {
        cout << "Не существует обратной матрицы." << endl;
    } else {
        // Вычисляем время выполнения
        t = clock() - t;
        double elapsed_time = t * 1.0 / CLOCKS_PER_SEC;
        cout << "Время работы алгоритма в секундах: " << elapsed_time << endl;

        // Сохраняем решение в файл
        ofstream outfile(output_file);
        for (int i = 0; i < n; ++i) {
            outfile << x[i] << " " << u[i] << endl;
        }
        outfile.close();
        cout << "Решение сохранено в файл " << output_file << endl;

        // Добавим запись в файл с таблицей времени выполнения
        ofstream timing_file("timing_results.txt", ios::app);
        timing_file << "h = " << h << ", Время выполнения: " << elapsed_time << " секунд" << endl;
        timing_file.close();
    }
    

    // Освобождаем выделенную память
    delete[] matrix;
    delete[] b;
    delete[] x;
    delete[] u;
    delete[] L;
    delete[] U;
    

    return 0;
}
