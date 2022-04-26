from random import randint
import math
import statistics
import numpy as np
import matplotlib.pyplot as plt
from MyUtils import *

def median_filter(data: list, N: int):
    filtered_data = []
    for i in range(1, len(data) - 1, 1):
        neighbours = data[i - math.floor(N/2) : i + math.ceil(N/2)]
        neighbours.sort()
        filtered_data.append(statistics.median(neighbours))
    return filtered_data

x = np.linspace(0, 50, 50)
randdata = [randint(0,100) for _ in range(50)]

plt.plot(x, randdata, c = "green")
plt.plot(x[1:-1], median_filter(randdata, 3), c = "red")
plt.plot(x, moving_average(randdata, 3), c = "blue")
plt.savefig("./src/filter_test_results.png")
plt.show()
