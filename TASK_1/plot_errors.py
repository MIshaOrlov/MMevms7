import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def calculate_u_true(x):
    #return np.sin(x * np.pi) - 2 * x*x*x + (3 + np.pi) * x*x - np.pi * x
    return - 2 * x*x*x + 3* x*x
    #return np.sin(x * 3.14159) - 2 * x*x*x + (3 + 3.14159) * x*x - 3.14159 * x

def calculate_errors(x, u, u_true,h):
    #print('беброчка')
    #print(u_true)
    #print(u)
    errors = u - u_true
    max_error = np.max(np.abs(errors))
    mse = np.sqrt(np.sum(errors**2)*h)
    return max_error, mse

def calculate_h(x):

    h = np.diff(x)
    return np.mean(h) if len(h) > 0 else 0  # Среднее значение h для файла

def process_files(folder_path):

    h_values = []
    max_errors = []
    mse_errors = []
    
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path) and filename.endswith('.txt'):
            data = np.loadtxt(file_path)  # Предполагаем, что формат файла простой
            x, u, u_true = data[:, 0], data[:, 1] , data[:, 2]
            # Вычисляем u_true для данного x
            #u_true = calculate_u_true(x)
            #print(x,u,u_true)
            #print('gjgf')
            h = calculate_h(x)
            max_error, mse_error = calculate_errors(x, u, u_true,h)
            
            h_values.append(h)
            max_errors.append(max_error)
            mse_errors.append(mse_error)

    # Сортировка данных по h в порядке убывания
    sorted_indices = np.argsort(h_values)[::-1]  # Индексы для сортировки по убыванию
    h_values = np.array(h_values)[sorted_indices]
    max_errors = np.array(max_errors)[sorted_indices]
    mse_errors = np.array(mse_errors)[sorted_indices]

    return h_values, max_errors, mse_errors

def calculate_convergence_order(errors, h_values):
    p_values = []
    for i in range(len(h_values) - 1):
        if errors[i] != 0 and errors[i + 1] != 0:  # Защита от деления на ноль
            p = np.log(errors[i] / errors[i + 1]) / np.log(h_values[i] / h_values[i + 1])
            p_values.append(p)
        else:
            p_values.append(np.nan)  # Если ошибка нулевая, порядок не определен
    p_values.append(np.nan)  # Нет значения порядка для последнего элемента
    return np.array(p_values)

def save_tables(h_values, max_errors, mse_errors, p_max, p_mse):
    N_values = 1 / h_values  # Вычисление N как 1/h

    # Создаем таблицу данных
    df = pd.DataFrame({
        'h': h_values,
        'N = 1/h': N_values,
        'Max Error': max_errors,
        'MSE Error': mse_errors,
        'Order of Convergence (Max Error)': p_max,
        'Order of Convergence (MSE Error)': p_mse
    })

    # Сохранение таблиц для каждого графика
    df[['h', 'Max Error', 'Order of Convergence (Max Error)']].to_csv('table_max_error_vs_h.txt', index=False, sep='\t')
    df[['h', 'MSE Error', 'Order of Convergence (MSE Error)']].to_csv('table_mse_error_vs_h.txt', index=False, sep='\t')
    df[['N = 1/h', 'Max Error', 'Order of Convergence (Max Error)']].to_csv('table_max_error_vs_N.txt', index=False, sep='\t')
    df[['N = 1/h', 'MSE Error', 'Order of Convergence (MSE Error)']].to_csv('table_mse_error_vs_N.txt', index=False, sep='\t')

def plot_errors(h_values, max_errors, mse_errors, p_max, p_mse):
    N_values = 1 / h_values  # Вычисление N как 1/h

    plt.figure(figsize=(12, 12))

    plt.subplot(2, 2, 1)
    plt.plot(h_values, max_errors, 'o-', label='Max Error')
    #plt.yscale('log')  # Логарифмическая шкала по оси h
    plt.xlabel('h')
    plt.ylabel('Max Error')
    plt.title('Зависимость максимальной ошибки от h')
    plt.legend()

    plt.subplot(2, 2, 2)
    plt.plot(h_values, mse_errors, 'o-', label='MSE Error', color='orange')
    plt.yscale('log')  # Логарифмическая шкала по оси h
    plt.xlabel('h')
    plt.ylabel('MSE Error')
    plt.title('Зависимость среднеквадратичной ошибки от h')
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(N_values, max_errors, 'o-', label='Max Error')
    plt.xlabel('N = 1/h')
    plt.ylabel('Max Error')
    plt.title('Зависимость максимальной ошибки от N')
    plt.legend()


    plt.subplot(2, 2, 4)
    plt.plot(N_values, mse_errors, 'o-', label='MSE Error', color='orange')
    plt.xlabel('N = 1/h')
    plt.ylabel('MSE Error')
    plt.title('Зависимость среднеквадратичной ошибки от N')
    plt.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <folder_path>")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    h_values, max_errors, mse_errors = process_files(folder_path)
    

    p_max = calculate_convergence_order(max_errors, h_values)
    p_mse = calculate_convergence_order(mse_errors, h_values)
    
    save_tables(h_values, max_errors, mse_errors, p_max, p_mse)  # Сохранение таблиц данных
    plot_errors(h_values, max_errors, mse_errors, p_max, p_mse)  # Построение графиков
