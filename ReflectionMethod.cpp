#include "ReflectionMethod.h"

int ReflectionInverse(int n, double* matrix, double* adjoint){
    double sum = 0.0; // s_k
    double aNorm = 0.0; // ||a_1ˆ(k-1)||
    double xNorm = 0.0; // ||xˆ(k)||
    double Uk = 0.0; // U(x)
    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            adjoint[i*n + j] = (i == j)*1.0;
        }
    }
    
    for(int k = 0; k < n; k++){
        
        sum = 0.0;
        
        //(12) Формула для подсчета s_k
        for(int j = k + 1; j < n ;  j++  ){
            sum += matrix[j*n + k] * matrix[j*n + k];
            
        }
        //(13) Формула подсчета нормы a
        aNorm = sqrt(matrix[k*n + k]*matrix[k*n + k] + sum);
        
        //Проверка случая 0 определителя
        if(aNorm < 1e-32 ){
            return -1;
        }
        
        matrix[k*n + k] -= aNorm;
        
        //(14) Формула подсчета нормы x
        xNorm = sqrt((matrix[k*n + k])*(matrix[k*n + k]) + sum);
        
        if (xNorm < 1e-100){
            matrix[k * n + k] += aNorm;
            continue;
        }
        
        //(16) Подсчет искомого вектора
        for(int j = k; j < n; j++){
            matrix[j*n + k] = 1.0 * matrix[j*n + k] / xNorm;
        }
        
      
        
        //(10) переход к новой матрице домножением на U(x) для матрицы исходной
        for(int i = k; i < n;i++){
            Uk = 0.0;
            for(int j = k; j < n; j++){
                Uk += matrix[j * n + k] * matrix[j * n + i];
            }
            for(int j = k; j < n; j++){
                matrix[j * n + i] -= 2.0 * Uk * matrix[j * n + k];
            }
        }
        
        //(10) переход к новой матрице домножением на U(x)  для матрицы доп
        for (int i = 0; i < n; i++ )
        {
            Uk = 0.0;
            for (int j = k; j < n; j++){
                Uk += matrix[j * n + k] * adjoint[j * n + i];
            }
            for (int j = k; j < n; j++){
                adjoint[j * n + i] -=  2.0 * Uk * matrix[j * n + k];
            }
                
        }
        
        matrix[k * n + k] = aNorm;

    }
    //Обратный ход Гаусса
    for (int i = 0; i < n; i++)
        for (int j = n - 1; j >= 0; j--){
            Uk = adjoint[j * n + i];
            for (int k = j + 1; k < n; k++){
                Uk -= matrix[j * n + k] * adjoint[k * n + i];
                if( Uk > 1e+13){
                    return -1;
                }
            }
            adjoint[j * n + i] = Uk / matrix[j * n + j];
            
        }
    
    return 1;
}


int multiplyMatrixByVector(int n, double* matrix, double* b, double* result) {
    // Инициализация результирующего вектора нулями
    for (int i = 0; i < n; ++i) {
        result[i] = 0.0;
    }

    // Перемножение матрицы на вектор
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += matrix[i * n + j] * b[j];
        }
    }
    
    return 1;
}


double ResidualCalc( double* a,  double* b, double* result, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += a[i * n + k] * b[k * n + j];
            }
            result[i * n + j] = sum;
            if( i == j) {
                result[i * n + j] -= 1;
                
            }
            std::cout<<result[i * n + j]<<" ";
           
        }
        std::cout<<std::endl;
        
    }
    return calculateRowNorm(result,n);
}



// Функция для LU-разложения пятидиагональной матрицы
int LU_decomposition(int N, const double* A, double* L, double* U) {
    // Инициализация L и U как нулевых матриц
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            L[i * N + j] = (i == j) ? 1.0 : 0.0;  // Единичная диагональ для L
            U[i * N + j] = 0.0;                   // Нулевая матрица для U
        }
    }

    // Выполнение LU-разложения
    for (int i = 0; i < N; ++i) {
        // Заполнение верхней треугольной матрицы U
        for (int j = i; j < std::min(i + 3, N); ++j) {  // Пятидиагональная структура
            U[i * N + j] = A[i * N + j];
            for (int k = 0; k < i; ++k) {
                U[i * N + j] -= L[i * N + k] * U[k * N + j];
            }
        }

        // Заполнение нижней треугольной матрицы L
        for (int j = i + 1; j < std::min(i + 3, N); ++j) {  // Пятидиагональная структура
            L[j * N + i] = A[j * N + i];
            for (int k = 0; k < i; ++k) {
                L[j * N + i] -= L[j * N + k] * U[k * N + i];
            }
            L[j * N + i] /= U[i * N + i];
        }
    }
    
    return 1;
}


int solve_LU(int N, const double* L, const double* U, double* X, const double* B) {
    double* Y = new double[N];  // Временный массив для решения LY = B

    // Прямой ход: решаем LY = B
    for (int i = 0; i < N; ++i) {
        Y[i] = B[i];
        for (int j = 0; j < i; ++j) {
            Y[i] -= L[i * N + j] * Y[j];
        }
        Y[i] /= L[i * N + i];  // Поскольку L имеет единичную диагональ, это просто Y[i] = Y[i]
    }

    // Обратный ход: решаем UX = Y
    for (int i = N - 1; i >= 0; --i) {
        X[i] = Y[i];
        for (int j = i + 1; j < N; ++j) {
            X[i] -= U[i * N + j] * X[j];
        }
        X[i] /= U[i * N + i];
    }

    delete[] Y;  // Освобождение памяти для временного массива
    return 1;
}
