import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def calculate_u_true(x):
    """Вычисляет истинные значения u_true по формуле."""
    return -2 * x**3 + 3 * x**2

def calculate_errors(x, u, u_true):
    """Рассчитывает максимальную и среднеквадратичную ошибки."""
    errors = u - u_true
    max_error = np.max(np.abs(errors))
    mse = np.sqrt(np.mean(errors**2))
    return max_error, mse

def calculate_h(x):
    """Рассчитывает шаг h как разность между соседними значениями x."""
    h = np.diff(x)
    return np.mean(h) if len(h) > 0 else 0  # Среднее значение h для файла

def process_files(folder_path):
    """Обрабатывает файлы и собирает ошибки и шаги."""
    h_values = []
    max_errors = []
    mse_errors = []
    
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path) and filename.endswith('.txt'):
            data = np.loadtxt(file_path)  # Предполагаем, что формат файла простой
            x, u = data[:, 0], data[:, 1]
            
            # Вычисляем u_true для данного x
            u_true = calculate_u_true(x)
            max_error, mse_error = calculate_errors(x, u, u_true)
            h = calculate_h(x)
            
            h_values.append(h)
            max_errors.append(max_error)
            mse_errors.append(mse_error)

    # Сортировка данных по h в порядке убывания
    sorted_indices = np.argsort(h_values)[::-1]  # Индексы для сортировки по убыванию
    h_values = np.array(h_values)[sorted_indices]
    max_errors = np.array(max_errors)[sorted_indices]
    mse_errors = np.array(mse_errors)[sorted_indices]

    return h_values, max_errors, mse_errors

def plot_errors(h_values, max_errors, mse_errors):
    """Строит графики зависимости ошибок от h и от N=1/h с сортировкой по убыванию h."""
    # Вычисление N как 1/h
    N_values = 1 / h_values

    plt.figure(figsize=(12, 12))

    # График зависимости максимальной ошибки от h
    plt.subplot(2, 2, 1)
    plt.plot(h_values, max_errors, 'o-', label='Max Error')
    plt.xlabel('h')
    plt.ylabel('Max Error')
    plt.title('Зависимость максимальной ошибки от h')
    plt.legend()

    # График зависимости среднеквадратичной ошибки от h
    plt.subplot(2, 2, 2)
    plt.plot(h_values, mse_errors, 'o-', label='MSE Error', color='orange')
    plt.xlabel('h')
    plt.ylabel('MSE Error')
    plt.title('Зависимость среднеквадратичной ошибки от h')
    plt.legend()

    # График зависимости максимальной ошибки от N = 1/h
    plt.subplot(2, 2, 3)
    plt.plot(N_values, max_errors, 'o-', label='Max Error')
    plt.xlabel('N = 1/h')
    plt.ylabel('Max Error')
    plt.title('Зависимость максимальной ошибки от N')
    plt.legend()

    # График зависимости среднеквадратичной ошибки от N = 1/h
    plt.subplot(2, 2, 4)
    plt.plot(N_values, mse_errors, 'o-', label='MSE Error', color='orange')
    plt.xlabel('N = 1/h')
    plt.ylabel('MSE Error')
    plt.title('Зависимость среднеквадратичной ошибки от N')
    plt.legend()

    plt.tight_layout()
    plt.show()

# Основной код для запуска из командной строки
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <folder_path>")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    h_values, max_errors, mse_errors = process_files(folder_path)
    plot_errors(h_values, max_errors, mse_errors)
