from scipy.signal import find_peaks, medfilt
import numpy as np
import statistics
import math

# from a matrix-like dataset, extracts the nth column
def extract_column(data:list, n):
    ext_data = []
    for i in range(len(data)):
        ext_data.append(data[i][n])
    return ext_data

# returns the mobile mean curve of the dataset given
def moving_average(data:list, N:int):
    offset = math.floor(N/2)
    averages = data[:offset]
    i = offset
    while i < len(data) - offset:
        window = data[i - offset : i + offset]
        window_average = sum(window) / N
        averages.append(window_average)
        i += 1
    averages.extend(data[len(data) - offset :])
    return averages

# applies median filter to dataset in given ranges (if not given, the gets applied to the whole domain)
# ranges follow this structure: [[a1,b1], [a2,b2], [a3,b3], ...]
def median_filter(data:list, N:int, ranges:list=None):
    if ranges is None:
        return medfilt(data, N).tolist()
    else:
        # [!] check if ranges intersect or are not in ascending order
        locally_filtered_data = data[:ranges[0][0]]
        for i in range(len(ranges)):
            locally_filtered_data.extend( medfilt(data[ranges[i][0]:ranges[i][1]], N).tolist() )
            if i + 1 < len(ranges):
                locally_filtered_data.extend( data[ranges[i][1]:ranges[i+1][0]] )
        locally_filtered_data.extend( data[ranges[-1][1]:] )
        if len(data) != len(locally_filtered_data):
            # print(data)
            # print(locally_filtered_data)
            raise Exception(f"Data corrupted: found lengths {len(data)} and {len(locally_filtered_data)}")
        return locally_filtered_data

# calculate the derivative function and returns it for the whole domain
def calc_derivative(data:list):
    derivative = []
    for i in range(len(data)):
        inext = i + 1
        if inext >= len(data):
            break
        else:
            # calculate incremental ratio between current data point and next (the last point in the set gets excluded):
            # d_dose / d_pos = (dose[i] - dose[inext]) / (pos[i] - pos[inext])
            derivative.append(
                (data[i][1] - data[inext][1]) / (data[i][0] - data[inext][0])
            )
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
def normalize_data(data:list):
    dose_max = max(extract_column(data, 1))
    norm_data = []
    for v in data:
        norm_data.append( (v / dose_max) * 100 )
    return norm_data
