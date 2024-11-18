# Компилятор и флаги
CXX = g++
CXXFLAGS = -Wall -std=c++11

# Файлы
SRC_DIR = .
BIN_DIR = ./bin
OUT_DIR = ./output
IMG_DIR = ./images

# Исходные файлы
SRCS = $(SRC_DIR)/main.cpp $(SRC_DIR)/RungeKutta.cpp
TARGET = $(BIN_DIR)/solver

# Значения h, для которых будем строить графики и замерять время
#HS = 0.01 0.005 0.0025 0.001 0.0005 0.00025 0.0001 0.00005 0.00001 0.000128 0.000256 0.000512 0.001024
#HS = 0.02 0.04 0.08 0.16 0.32 0.64
HS = 0.1 0.05 0.025 0.01 0.005 0.0025 0.001 
# Файл для записи времени выполнения
TIMING_FILE = timing_results.txt

# Цель по умолчанию
all: run_all_plots

# Компиляция программы
$(TARGET): $(SRCS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

# Запуск программы с разными значениями h, построение графиков и запись времени
run_all_plots: $(TARGET)
	@mkdir -p $(OUT_DIR)
	@mkdir -p $(IMG_DIR)
	@echo "h и время выполнения" > $(TIMING_FILE)
	@for h in $(HS); do \
		echo "Запуск с шагом $$h и параметром 2..."; \
		$(TARGET) $$h 2 $(OUT_DIR)/solution_$$h.txt; \
		python3 plot_solution.py $(OUT_DIR)/solution_$$h.txt $(IMG_DIR)/plot_$$h.png $$h; \
	done
	@echo "Графики построены и сохранены в $(IMG_DIR)"
	python3 plot_timing.py $(TIMING_FILE) $(IMG_DIR)/timing_plot.png
	python3 plot_errors.py $(OUT_DIR)




# Запуск программы в тестовом режиме (параметр 1) с фиксированным h = 0.0001
run_test: $(TARGET)
	@echo "Запуск программы в тестовом режиме с h = 0.0001 и параметром 1..."
	$(TARGET) 0.0001 1 $(OUT_DIR)/test_solution.txt
	python3 plot_solution.py $(OUT_DIR)/test_solution.txt $(IMG_DIR)/test_plot.png 0.0001
	@echo "График для тестового режима сохранен в $(IMG_DIR)/test_plot.png"
	

# Очистка всех сгенерированных файлов
clean:
	rm -rf $(BIN_DIR) $(OUT_DIR) $(IMG_DIR) $(TIMING_FILE)
