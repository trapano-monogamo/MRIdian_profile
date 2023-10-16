from scipy.signal import medfilt
import numpy as np
import math


# from a matrix-like dataset, extracts the nth column
def extract_column(data: list, n):
    ext_data = []
    for i in range(len(data)):
        ext_data.append(data[i][n])
    return ext_data


def chi_squared(observed: list, expected: list) -> float:
    acc = 0
    for k in range(1, len(observed)):
        acc += (observed[k] - expected[k]) ** 2 / expected[k]
    return acc


def windowed_chi_squared(observed: list, expected: list, a: int, b: int) -> float:
    acc = 0
    for k in range(a, b):
        acc += (observed[k] - expected[k]) ** 2 / abs(expected[k])
    return acc


def windowed_reduced_chi_squared(observed: list, expected: list, a: int, b: int) -> float:
    acc = 0
    _acc = 0
    for k in range(a, b):
        acc += (observed[k] - expected[k]) ** 2 / abs(expected[k])
        _acc += (observed[k] - expected[k]) ** 2
    print(f"N = {abs(b-a)}\no-e = {(observed[(a - b) // 3] - expected[(a-b) // 3]) ** 2}\nacc = {acc}\nchi2 = {acc / (abs(b-a) - 1)}\n_acc = {_acc / abs(b-a)}")
    return acc / (abs(b - a) - 1)


# returns the mobile mean curve of the dataset given
def moving_average(data: list, N: int):
    offset = math.floor(N / 2)
    averages = data[0:offset]
    window_average = 0
    i = offset
    while i < len(data) - offset:
        # take into account current point, which is at (i - offset) + 1  =>  therefore last point is i + offset + 1
        window = data[i - offset:i + offset + 1]
        window_average = sum(window) / N
        averages.append(window_average)

        # print(f"{i}: {window}, {window_average}")

        if i == offset:  # first iteration
            # fill head of dataset with first modified value
            for n in range(0, offset, 1):  # 0 -> offset (1)
                averages[n] = window_average

        if i == len(data) - 1 - offset:  # last iteration
            # fill tail with last modified value
            averages.extend([window_average for _ in range(offset)])
        i += 1
    return averages


# applies median filter to dataset in given ranges (if not given, the gets applied to the whole domain)
# ranges follow this structure: [[a1,b1], [a2,b2], [a3,b3], ...]
def median_filter(data: list, N: int, ranges: list = None):
    if ranges is None:
        return medfilt(data, N).tolist()
    else:
        # [!] check if ranges intersect or are not in ascending order
        locally_filtered_data = data[:ranges[0][0]]
        for i in range(len(ranges)):
            locally_filtered_data.extend(medfilt(data[ranges[i][0]:ranges[i][1]], N).tolist())
            if i + 1 < len(ranges):
                locally_filtered_data.extend(data[ranges[i][1]:ranges[i + 1][0]])
        locally_filtered_data.extend(data[ranges[-1][1]:])
        if len(data) != len(locally_filtered_data):
            # print(data)
            # print(locally_filtered_data)
            raise Exception(f"Something went wrong: found lengths {len(data)} and {len(locally_filtered_data)}")
        return locally_filtered_data


# calculate the derivative function and returns it for the whole domain
def calc_derivative(x: list, y: list):
    derivative = []
    for i in range(len(x)):
        inext = i + 1
        if inext >= len(x):
            break
        else:
            # calculate incremental ratio between current data point and next (the last point in the set gets excluded):
            # dy / dx = (dose[i] - dose[inext]) / (pos[i] - pos[inext])
            derivative.append(
                (y[i] - y[inext]) / (x[i] - x[inext])
            )
    derivative.append(derivative[-1])
    return derivative


def find_intersections(f: list, g: list):
    if len(f) != len(g):
        raise Exception(f"Lists must have the same length (found {len(f)} and {len(g)})")
    else:
        intersections = []
        x = 0
        diff_sign = np.sign(0)
        first_iteration = True
        while x < len(f):
            new_diff_sign = np.sign(f[x] - g[x])
            if not first_iteration:
                if new_diff_sign.astype(int) != diff_sign.astype(int):
                    intersections.append(x)
            first_iteration = False
            diff_sign = new_diff_sign
            x += 1
        return intersections


# normalizes dataset
def normalize_data(data: list):
    dose_max = max(extract_column(data, 1))
    norm_data = []
    for v in data:
        norm_data.append((v / dose_max) * 100)
    return norm_data


def transpose_table(table):
    result = []
    for i in range(len(table[0])):
        row = []
        for item in table:
            row.append(item[i])
        result.append(row)
    return result


def gauss(x, amp, cen, wid):
    return amp * np.exp(-(x - cen) ** 2.0 / (2 * wid ** 2))


def gauss_first_derivative(x, amp, cen, wid):
    return -amp * (x - cen) * np.exp(-(x - cen) ** 2.0 / (2 * wid ** 2.0)) / wid ** 2.0


def gauss_second_derivative(x, amp, cen, wid):
    return amp * (x ** 2.0 - 2 * cen * x - wid ** 2.0 + cen ** 2.0) * np.exp(-(x - cen) ** 2.0 / (2 * wid ** 2.0)) / wid ** 4.0


def lerp(v0, v1, t):
    # return v0 + t * (v1 - v0)
    return (1 - t) * v0 + t * v1


def max_and_min_in_range(_data: list, _begin: int = None, _end: int = None):
    begin = _begin if _begin else 0
    end = _end if _end else len(_data)
    dmax = max(_data[begin:end])
    dmin = min(_data[begin:end])
    dmaxi = _data.index(dmax, begin, end)
    dmini = _data.index(dmin, begin, end)
    return (dmax, dmaxi, dmin, dmini)


def continuous_max_and_min_in_range(self, f, params: list, n: int, begin: float = 0.0, end: float = 1.0):
    domain = end - begin
    dmax = f(begin, *params)
    dmaxi = 0.0
    dmin = f(begin, *params)
    dmini = 0.0
    # iterate through n points from _begin to _end
    x = begin
    while x < end:
        f_x = f(x, *params)
        # if new max, update and save position
        if dmax < f_x:
            dmax = f_x
            dmaxi = x
        # if new min, update and save position
        if dmin > f_x:
            dmin = f_x
            dmini = x
        x += (domain / n)
    return (dmax, dmaxi, dmin, dmini)
