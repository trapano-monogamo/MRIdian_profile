# TODO: remove \t from Profile->raw_data, and avoid those horrible "\t\tSOMETHING"...

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, medfilt
from MyUtils import *


log = False




# --- Scan ---
# stores scan data and parameters in a dictionary.

class Scan:
   scan_num: int
   fields: dict
   data: list

   derivative: list
   filtered_derivative: list
   inflection_points: list

   profile_out_dir: str

   def __init__(self, _scan_num, _raw_data, _begin, _end, _profile_out_dir):
      self.scan_num = _scan_num
      self.fields = {}
      self.data = []
      self.derivative = []
      self.filtered_derivative = []
      self.inflection_points = []
      self.profile_out_dir = _profile_out_dir

      # log
      if log:
         print(f"Reading profile {self.scan_num}:")

      # loop through the region, as long as "BEGIN_DATA" isn't encountered, keep collecting fields,
      # once it is encountered, collect floats
      i = _begin + 1
      while i != _end:
         # if it's a field, separate field name and value (they're separated with an equal sign), and put them in the fields dict
         if _raw_data[i].find("BEGIN_DATA") == -1:
            tokens = _raw_data[i].replace("\t", "").split("=")
            self.fields[tokens[0]] = tokens[1]
         # otherwise collect the numbers in the current line as pairs and put them in the data list
         else:
            j = i + 1
            while j != _end - 1:
               # removing empty strings (once separators) from list and cast the rest to floats
               nums = [e for e in _raw_data[j].split("\t") if e]
               self.data.append([float(nums[0]), float(nums[1])])
               j += 1
            # stop parsing
            break
         i += 1

      self.output_plot()


   def collect_fields(self):
      pass


   def find_inflection_points(self):
      if log:
         print(f"\tInflection Points...", end="")

      # filtered_data = medfilt(self.extract_data_column(1), 3)
      # for i in range(len(filtered_data)):
      #    self.data[i][1] = filtered_data[i]

      self.derivative = calc_derivative(self.data)

      data_average = [np.mean(extract_column(self.data, 1)) / 3.0 for _ in range(len(self.data))]
      # data_average = moving_average(extract_column(self.data, 1), 50)
      intersections = find_intersections(extract_column(self.data, 1), data_average)
      # self.derivative = medfilt(self.derivative, 5).tolist()
      self.derivative = median_filter(self.derivative, 5, [[0,intersections[0]], [intersections[1], len(self.derivative)]])

      pmax = max(self.derivative)
      pmaxi = self.derivative.index(pmax)
      pmin = min(self.derivative)
      pmini = self.derivative.index(pmin)

      print(f"{self.scan_num}: {pmaxi > intersections[0]}")

      self.inflection_points = [
         [pmaxi, pmax], # positive peak
         [pmini, pmin],  # negative peak
      ]

      if log:
         dose_data = extract_column(self.data, 1)
         print(f"\t\tDone", end="\n")
         print(f"\t- {pmaxi} ; {pmini}", end="\n")
         print(f"\t- ({pmax:.3f} : {dose_data[pmaxi]:.3f}) ; ({pmin:.3f} : {dose_data[pmini]:.3f})", end="\n")
      

   def output_plot(self):
      self.find_inflection_points()

      if log:
         print(f"\tFinal Plot...", end="")

      x = np.linspace(0, len(self.derivative), len(self.derivative))
      y = self.derivative
      # y2 = moving_average(extract_column(self.data[:-1], 1), 50)
      y2 = [np.mean(self.data) / 3.0 for _ in range(len(self.data) - 1)]

      dose_data = extract_column(self.data[:-1], 1)

      fig, ax = plt.subplots()
      ax.plot(x, dose_data, c = "blue")
      ax.plot(x, y, c = "red")
      ax.plot(x, y2, c = "purple")
      ax.scatter(
         [self.inflection_points[0][0], self.inflection_points[1][0]],
         [self.inflection_points[0][1], self.inflection_points[1][1]],
         c = "green")
      plt.savefig(f"{self.profile_out_dir}/{self.scan_num}.png")
      plt.close(fig)

      if log:
         print(f"\t\t\tDone", end="\n")




# --- Profile ---
# reads a file and extracts scan data present in that file
# by parsing scan by scan each field and collecting numerical data,
# organised in DataScans instances.

class Profile:
   name: str
   scans: list

   def __init__(self, file_path: str, out_dir: str):
      self.scans = []
      self.name = file_path.split("/")[-1].replace(".mcc", "")
      raw_data = []

      # log
      if log:
         print(f"O===== Reading [{self.name}] =====O")
      print(f"{self.name}")

      # create output profile subdirectory
      profile_out_dir = f"{out_dir}/{self.name}/"
      if not os.path.exists(profile_out_dir):
         os.mkdir(profile_out_dir)

      # read file
      with open(file_path, "r") as f:
         raw_data = f.read().split("\n")
         
      # parsing
      i = 0
      end_region_index = 0
      while i < len(raw_data):
         # when a BEGIN_SCAN is found, highlight its region up to END_SCAN, and build a DataScan out of it
         if "\tBEGIN_SCAN" in raw_data[i].split(' '):
            scan_num = int(raw_data[i].split(" ")[-1])

            # OMG!!! what are we, babies?! c'mon...
            if scan_num < 10:
               end_region_index = raw_data.index(f"\tEND_SCAN  {scan_num}")
            else:
               end_region_index = raw_data.index(f"\tEND_SCAN {scan_num}")

            self.scans.append(Scan(scan_num, raw_data, i, end_region_index, profile_out_dir))

            # skip content in between
            i = end_region_index
         i += 1

      print("")



# --- Cacher ---
# stores all Profiles present in a folder

class Cacher:
   profiles: list

   def __init__(self, res_dir: str, out_dir: str):
      self.profiles = []

      # create filelist
      filelist = []
      for root, dirs, files in os.walk(res_dir):
         for _file in files:
            current_file_path = os.path.join(root,_file).replace("\\", "/")
            filelist.append(current_file_path)

      # create profile out of each file in filelist
      for f in filelist:
         self.profiles.append(Profile(f, out_dir))
