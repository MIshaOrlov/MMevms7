#include "RungeKutta.h"




double none(double x){
    return 0;
}




// Функция правой части системы ОДУ для неоднородного уравнения
void ODESystemNonHomogeneous(double x, const double y[], double dydx[]) {
    // y[0] = u
    // y[1] = u'
    // y[2] = u''
    // y[3] = u'''
    dydx[0] = y[1];
    dydx[1] = y[2];
    dydx[2] = y[3];
    dydx[3] = std::exp(x * x) - y[3] - std::sin(x) * y[0];
}



// Метод Рунге-Кутты второго порядка для системы ОДУ
void RungeKutta2( double (*func)(double), double h, const double y0[], double x_vals[], double y_vals[], int N_points, int dim) {

    // Инициализация массивов
    double *y = new double[dim];
    double *k1 = new double[dim];
    double *k2 = new double[dim];

    x_vals[0] = 0;
    for (int j = 0; j < dim; ++j) {
        y_vals[0 * dim + j] = y0[j];
    }

    for (int i = 0; i < N_points - 1; ++i) {
        
        for (int j = 0; j < dim; ++j) {
            y[j] = y_vals[i * dim + j];
        }

        //ODE(x_vals[i], y, k1);
        k1[0] = y[1];
        k1[1] = y[2];
        k1[2] = y[3];
        k1[3] = func(x_vals[i]) - y[3] - std::sin(x_vals[i]) * y[0];

    
        for (int j = 0; j < dim; ++j) {
            y[j] += h * k1[j];
        }
        //ODE(x_vals[i] + h, y, k2);
        k2[0] = y[1];
        k2[1] = y[2];
        k2[2] = y[3];
        k2[3] = func(x_vals[i]) - y[3] - std::sin(x_vals[i]) * y[0];

        for (int j = 0; j < dim; ++j) {
            y_vals[(i + 1) * dim + j] = y_vals[i * dim + j] + h * 0.5 * (k1[j] + k2[j]);
            y[j] = y_vals[i * dim + j];
        }
    }

    // Освобождение памяти
    delete[] y;
    delete[] k1;
    delete[] k2;
}



// Метод Рунге-Кутты четвертого порядка для системы ОДУ

// Метод Рунге-Кутты четвертого порядка для системы ОДУ
void RungeKutta4(double (*func)(double), double h, const double y0[], double x_vals[], double y_vals[], int N_points, int dim) {

    // Инициализация массивов
    double *y = new double[dim];
    double *k1 = new double[dim];
    double *k2 = new double[dim];
    double *k3 = new double[dim];
    double *k4 = new double[dim];

    // Начальное значение x
    x_vals[0] = 0.0;
    // Установка начальных условий для y
    for (int j = 0; j < dim; ++j) {
        y_vals[0 * dim + j] = y0[j];
    }

    // Основной цикл для всех точек
    for (int i = 0; i < N_points - 1; ++i) {
        
        // Получение текущего значения y
        for (int j = 0; j < dim; ++j) {
            y[j] = y_vals[i * dim + j];
        }
        // Вычисление k1
        k1[0] = h * y[1];
        k1[1] = h * y[2];
        k1[2] = h * y[3];
        k1[3] = h * (func(x_vals[i]) - y[3] - std::sin(x_vals[i]) * y[0]);

        // Вычисление k2, используя значения x и y, сдвинутые на 0.5 * h
        double x_half = x_vals[i] + 0.5 * h;
        for (int j = 0; j < dim; ++j) {
            y[j] = y_vals[i * dim + j] + 0.5 * k1[j];
        }
        k2[0] = h * y[1];
        k2[1] = h * y[2];
        k2[2] = h * y[3];
        k2[3] = h * (func(x_half) - y[3] - std::sin(x_half) * y[0]);

        // Вычисление k3, используя значения x и y, сдвинутые на 0.5 * h
        for (int j = 0; j < dim; ++j) {
            y[j] = y_vals[i * dim + j] + 0.5 * k2[j];
        }
        k3[0] = h * y[1];
        k3[1] = h * y[2];
        k3[2] = h * y[3];
        k3[3] = h * (func(x_half) - y[3] - std::sin(x_half) * y[0]);

        // Вычисление k4, используя значения x и y, сдвинутые на полный шаг h
        double x_full = x_vals[i] + h;
        for (int j = 0; j < dim; ++j) {
            y[j] = y_vals[i * dim + j] + k3[j];
        }
        k4[0] = h * y[1];
        k4[1] = h * y[2];
        k4[2] = h * y[3];
        k4[3] = h * (func(x_full) - y[3] - std::sin(x_full) * y[0]);

        // Обновление значений y_vals с учетом взвешенной суммы коэффициентов
        for (int j = 0; j < dim; ++j) {
            y_vals[(i + 1) * dim + j] = y_vals[i * dim + j] + (1.0 / 6.0) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
        }

        // Обновляем x на следующий шаг
        x_vals[i + 1] = x_vals[i] + h;
    }

    // Освобождение памяти
    delete[] y;
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
}


int ShooterMethod(double (*func)(double),int N_points,int dim,double h, double* y_vals, double* x_vals, double* u_vals,double* phi1_full,double* phi2_full,double* psi_full){
    // Нарезка сетки по x
    for(int i =0; i< N_points; i++){
        x_vals[i] = i * h;
    }
    
    // Решение для phi1 (фундаментальное решение 1)
    double y0_phi1[4] = {0.0, 0.0, 1.0, 0.0}; // u(0)=0, u'(0)=0, u''(0)=1, u'''(0)=0
    RungeKutta4(none, h, y0_phi1, x_vals, y_vals, N_points, dim);
    

    // Сохраняем значения phi1 и phi1'
    double phi1 = y_vals[(N_points - 1) * dim + 0];
    double phi1_prime = y_vals[(N_points - 1) * dim + 1];
    
    // Сохраняем решение phi1 для дальнейшего использования
    for (int i = 0; i < N_points; ++i) {
        for (int j = 0; j < dim; ++j) {
            phi1_full[i * dim + j] = y_vals[i * dim + j];
        }
    }

    // Решение для phi2 (фундаментальное решение 2)
    double y0_phi2[4] = {0.0, 0.0, 0.0, 1.0}; // u(0)=0, u'(0)=0, u''(0)=0, u'''(0)=1
    RungeKutta4(none, h, y0_phi2, x_vals, y_vals, N_points, dim);

    // Сохраняем значения phi2 и phi2'
    double phi2 = y_vals[(N_points - 1) * dim + 0];      // y[N][0]
    double phi2_prime = y_vals[(N_points - 1) * dim + 1]; // y[N][1]

    // Сохраняем решение phi2 для дальнейшего использования
    for (int i = 0; i < N_points; ++i) {
        for (int j = 0; j < dim; ++j) {
            phi2_full[i * dim + j] = y_vals[i * dim + j];
        }
    }

    // Решение для psi (частное решение)
    double y0_psi[4] = {0.0, 0.0, 0.0, 0.0}; // u(0)=0, u'(0)=0, u''(0)=0, u'''(0)=0
    RungeKutta4(func, h, y0_psi, x_vals, y_vals, N_points, dim);

    // Сохраняем значения psi и psi'
    double psi = y_vals[(N_points - 1) * dim + 0];      // y[N][0]
    double psi_prime = y_vals[(N_points - 1) * dim + 1]; // y[N][1]

    // Сохраняем решение psi для дальнейшего использования

    for (int i = 0; i < N_points; ++i) {
        for (int j = 0; j < dim; ++j) {
            psi_full[i * dim + j] = y_vals[i * dim + j];
        }
    }
    
    

    // Составляем систему линейных уравнений для нахождения C1 и C2
    // C1 * phi1(1) + C2 * phi2(1) + psi(1) = 1
    // C1 * phi1'(1) + C2 * phi2'(1) + psi'(1) = 0

    double b1 = 1.0 - psi;
    double b2 = 0.0 - psi_prime;

    // Матрица коэффициентов
    double a11 = phi1;
    double a12 = phi2;
    double a21 = phi1_prime;
    double a22 = phi2_prime;

    // Вычисляем определитель матрицы
    double det = a11 * a22 - a12 * a21;

    if (std::abs(det) < 1e-16) {
        return -1;
    }

    // Решаем систему методом Крамера
    double C1 = (b1 * a22 - b2 * a12) / det;
    double C2 = (a11 * b2 - a21 * b1) / det;

    
    for (int i = 0; i < N_points; ++i) {
        u_vals[i] = C1 * phi1_full[i * dim + 0] + C2 * phi2_full[i * dim + 0] + psi_full[i * dim + 0];
    }


    
    return 1;
}
