import matplotlib.pyplot as plt
import numpy as np

sequences = ['100', '1 000', '10 000', '100 000', '300 000']

time_single = [0.056,  0.0853, 0.943, 6.81, 16.24]
time_threads = [0.0189, 0.111, 1.072, 8.995, 21.28]
time_processes = [1.53, 1.49, 2.23, 3.81, 8.44]

x = np.arange(len(sequences))
width = 0.25

fig, ax = plt.subplots(figsize=(10, 6))

ax.bar(x - width, time_single, width, label='Однопоточный', color="#b80588")
ax.bar(x, time_threads, width, label='Многопоточный (4)', color="#3617bd")
ax.bar(x + width, time_processes, width, label='Многопроцессный (4)', color="#ab1010")

ax.set_title('Зависимость времени от количества последовательностей')
ax.set_xlabel('Количество последовательностей')
ax.set_ylabel('Время выполнения (сек)')
ax.set_xticks(x)
ax.set_xticklabels(sequences)
ax.legend()
ax.grid(axis='y', linestyle='--', alpha=0.7)

plt.show()


workers = ['1', '2', '4', '8', '12']

time_threads = [7.8, 9.04, 9.5, 10.16, 10.41]
time_processes = [12.96, 5.46, 3.92, 4.77, 5.42]

x = np.arange(len(workers))
width = 0.35

fig, ax = plt.subplots(figsize=(10, 6))

ax.bar(x - width/2, time_threads, width, label='Многопоточность', color="#31c713")
ax.bar(x + width/2, time_processes, width, label='Многопроцессность', color="#d1cd0e")

ax.set_title('Зависимость времени от числа потоков/процессов (файл 100k)')
ax.set_xlabel('Количество потоков / процессов')
ax.set_ylabel('Время выполнения (сек)')
ax.set_xticks(x)
ax.set_xticklabels(workers)
ax.legend()
ax.grid(axis='y', linestyle='--', alpha=0.7)

plt.show()
