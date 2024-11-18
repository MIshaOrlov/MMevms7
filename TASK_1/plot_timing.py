import matplotlib.pyplot as plt
import sys

def read_timing_data(filename):
    h_values = []
    times = []
    with open(filename, 'r') as file:
        next(file)  # Пропустить первую строку заголовка
        for line in file:
            # Извлекаем значение h и время выполнения из строки
            parts = line.split(',')
            h = float(parts[0].split('=')[1].strip())
            time = float(parts[1].split(':')[1].strip().split()[0])
            h_values.append(h)
            times.append(time)
    return h_values, times

def plot_timing(h_values, times, output_image):
    plt.plot(h_values, times, marker='o', linestyle='-', color='blue', label='Время выполнения')
    plt.xscale('log')
    plt.yscale('log')  
    plt.xlabel('Шаг h')
    plt.ylabel('Время выполнения (с)')
    plt.title('Зависимость времени выполнения от шага h')
    plt.grid(True)
    plt.legend()
    
    # Сохранение графика в файл
    plt.savefig(output_image)
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plot_timing.py <input_file> <output_image>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_image = sys.argv[2]

    h_values, times = read_timing_data(input_file)
    plot_timing(h_values, times, output_image)
