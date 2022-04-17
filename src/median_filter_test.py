from random import randint
import math
import statistics
import numpy as np
import matplotlib.pyplot as plt


# def median_filter(self, data: list, N: int):
#     offset = math.floor(N/2)
# 
#     filtered_data = []
#     filtered_data.append(data[:offset])
# 
#     for i in range(offset, len(data) - offset, 1):
#         neighbours = data[i - offset : i + offset]
#         neighbours.sort()
#         filtered_data.append(statistics.median(neighbours))
# 
#     filtered_data.append(data[-offset:])
# 
#     return filtered_data

def median_filter(data: list, N: int):
    filtered_data = []

    for i in range(1, len(data) - 1, 1):
        neighbours = data[i - math.floor(N/2) : i + math.ceil(N/2)]
        neighbours.sort()
        filtered_data.append(statistics.median(neighbours))

    return filtered_data

x = np.linspace(0, 50, 50)
data = [randint(0,100) for _ in range(50)]
filtered_data = median_filter(data, 3)

plt.plot(x, data, c = "green")
plt.plot(x[1:-1], filtered_data, c = "red")
plt.savefig("./src/median_filter_test_results.png")
plt.show()
