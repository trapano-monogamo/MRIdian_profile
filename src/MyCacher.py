# TODO: remove \t from Profile->raw_data, and avoid those horrible "\t\tSOMETHING"...

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from MyUtils import *




class TestSettings:
   test_preset: int
   test_method: str

   def __init__(self, tp , tm):
      self.test_preset = tp
      self.test_method = tm


# --- Scan ---
# stores scan data and parameters in a dictionary.

class Scan:
   # scan properties
   scan_num: int
   fields: dict
   pos_data: list
   dose_data: list

   # calculated data
   derivative: list
   second_derivative: list
   inflection_points: list

   # useful data
   profile_out_dir: str

   # test driven stuff babyyy
   test_settings: TestSettings
   test_presets: list

   def __init__(self, _scan_num, _raw_data, _begin, _end, _profile_out_dir, _test_settings: TestSettings):
      self.scan_num = _scan_num
      self.fields = {}
      self.pos_data = []
      self.dose_data = []
      self.derivative = []
      self.second_derivative = []
      self.inflection_points = []
      self.profile_out_dir = _profile_out_dir

      self.test_presets = [
         self.not_filtered_processing,
         self.dose_filtered_processing,
         self.first_derivative_filtered_processing,
         self.dose_first_derivative_filtered_processing,
         self.second_derivative_filtered_processing,
         self.both_derivatives_filtered_processing,
         self.all_filtered_processing
      ]
      self.test_settings = _test_settings

      # loop through the region, as long as "BEGIN_DATA" isn't encountered, keep collecting fields,
      # once it is encountered, collect floats
      temp_complete_data = []
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
               self.pos_data.append(float(nums[0]))
               self.dose_data.append(float(nums[1]))
               j += 1
            # stop parsing
            break
         i += 1

      # self.check_symmetry()
      # exit()

      self.produce_results()

   def check_symmetry(self):
      at_zero = self.pos_data.index(0.0)
      for i in range(0, len(self.pos_data)):
         if self.dose_data[i] != self.dose_data[(len(self.dose_data) - 1) - i]:
            print(f"pos: {self.pos_data[i]} - dose: {self.dose_data[i]} : {self.dose_data[-i]}")

   def apply_filter(self, _data:list):
      if self.test_settings.test_method == "moving_average":
         return moving_average(_data, 3)
      elif self.test_settings.test_method == "median_filter":
         # data_average = [np.mean(self.dose_data) / 3.0 for _ in range(len(self.dose_data))]
         # intersections = find_intersections(self.dose_data, data_average)
         return median_filter(_data, 3)
      elif self.test_settings.test_method == "spline":
         tck = splrep(self.pos_data[:-1], _data)
         return splev(self.pos_data[:-1], tck).tolist()
      else:
         raise Exception("Lol not yet")


   def not_filtered_processing(self):
      self.derivative = calc_derivative(self.pos_data, self.dose_data)
      self.second_derivative = calc_derivative(self.pos_data, self.derivative)
      

   def dose_filtered_processing(self):
      self.dose_data_raw = self.dose_data[:]
      self.dose_data = self.apply_filter(self.dose_data)
      self.derivative = calc_derivative(self.pos_data, self.dose_data)
      self.second_derivative = calc_derivative(self.pos_data, self.derivative)
      ranges = [[0,10], [len(self.second_derivative) - 10, len(self.second_derivative)]]
      self.second_derivative = median_filter(self.second_derivative, 3, ranges)


   def first_derivative_filtered_processing(self):
      self.derivative = calc_derivative(self.pos_data, self.dose_data)
      self.derivative = self.apply_filter(self.derivative)
      self.second_derivative = calc_derivative(self.pos_data, self.derivative)


   def dose_first_derivative_filtered_processing(self):
      self.dose_data = self.apply_filter(self.dose_data)
      self.derivative = calc_derivative(self.pos_data, self.dose_data)
      self.derivative = self.apply_filter(self.derivative)
      self.second_derivative = calc_derivative(self.pos_data, self.derivative)


   def second_derivative_filtered_processing(self):
      self.derivative = calc_derivative(self.pos_data, self.dose_data)
      self.second_derivative = calc_derivative(self.pos_data, self.derivative)
      self.second_derivative = self.apply_filter(self.second_derivative)


   def both_derivatives_filtered_processing(self):
      self.derivative = calc_derivative(self.pos_data, self.dose_data)
      self.derivative = self.apply_filter(self.derivative)
      self.second_derivative = calc_derivative(self.pos_data, self.derivative)
      self.second_derivative = self.apply_filter(self.second_derivative)



   def all_filtered_processing(self):
      self.dose_data = self.apply_filter(self.dose_data)
      self.derivative = calc_derivative(self.pos_data, self.dose_data)
      self.derivative = self.apply_filter(self.derivative)
      self.second_derivative = calc_derivative(self.pos_data, self.derivative)
      self.second_derivative = self.apply_filter(self.second_derivative)


   # Calculate derivative, than apply median filter to portions of the derivative to
   # eliminate peaks generated by experimental errors jumps in the original data.
   # These portions are identified as external intervals of the intersections between
   # the original data curve and a third of the mean (arbitrary value... ugh).
   # Once filtered, max and min of the derivative will be the functions high and low,
   # which are the derivative values at the inflection points positions.
   # To find the positions take the index of the peak in the derivative, and use it as index for
   # both derivative function and original data
   def find_inflection_points(self):

      # calculate and locally filter the derivative
      # self.dose_data = median_filter(self.dose_data, 5)
      # self.derivative = calc_derivative(self.pos_data, self.dose_data)
      # self.derivative = median_filter(self.derivative, 5, [[0,intersections[0]], [intersections[1], len(self.derivative)]])
      # self.derivative = median_filter(self.derivative, 5)
      # self.second_derivative = calc_derivative(self.pos_data[:-1], self.derivative)
      # self.second_derivative = median_filter(self.second_derivative, 5)
      self.test_presets[self.test_settings.test_preset]()

      # find peaks/valleys and their positions
      pmax = max(self.derivative)
      pmaxi = self.derivative.index(pmax)
      pmin = min(self.derivative)
      pmini = self.derivative.index(pmin)

      spmax = max(self.second_derivative[:len(self.second_derivative) // 2])
      spmaxi = self.second_derivative.index(spmax)
      spmin = min(self.second_derivative[:len(self.second_derivative) // 2])
      spmini = self.second_derivative.index(spmin)
      spmax2 = max(self.second_derivative[len(self.second_derivative) // 2:])
      spmaxi2 = self.second_derivative.index(spmax)
      spmin2 = min(self.second_derivative[len(self.second_derivative) // 2:])
      spmini2 = self.second_derivative.index(spmin)

      # save inflection points data: [pos, [derivative_value, data_value]]
      self.inflection_points = [
         [ (self.pos_data[pmaxi] + self.pos_data[pmaxi + 1]) / 2.0, [pmax, (self.dose_data[pmaxi] + self.dose_data[pmaxi+1]) / 2.0] ],
         [ (self.pos_data[pmini] + self.pos_data[pmini + 1]) / 2.0, [pmin, (self.dose_data[pmini] + self.dose_data[pmini+1]) / 2.0] ],
         [ self.pos_data[spmaxi], [spmax, self.dose_data[spmaxi]] ],
         [ self.pos_data[spmini], [spmin, self.dose_data[spmini]] ],
         [ self.pos_data[spmaxi2], [spmax2, self.dose_data[spmaxi2]] ],
         [ self.pos_data[spmini2], [spmin2, self.dose_data[spmini2]] ],
      ]

      # self.inflection_points = [
      #    # first derivative
      #    [
      #       (self.pos_data[pmaxi] + self.pos_data[pmaxi + 1]) / 2.0,
      #       [pmax, (self.dose_data[pmaxi] + self.dose_data[pmaxi + 1]) / 2.0]
      #    ],
      #    [
      #       (self.pos_data[pmini] + self.pos_data[pmini + 1]) / 2.0,
      #       [pmin, (self.dose_data[pmini] + self.dose_data[pmini + 1]) / 2.0]
      #    ],
      #    # second derivative (only on first half)
      #    [
      #       self.pos_data[spmaxi], [spmax, self.dose_data[spmaxi]]
      #    ],
      #    [
      #       self.pos_data[spmini], [spmin, self.dose_data[spmini]]
      #    ]
      # ]

      # print(f"pos: {self.inflection_points[0][0] == -self.inflection_points[1][0]}, deriv: {self.inflection_points[0][1][0] == -self.inflection_points[1][1][0]}, dose: {self.inflection_points[0][1][1] == self.inflection_points[1][1][1]}, filt: {pmaxi > intersections[0]},{pmini < intersections[1]}")

      # for i in range(len(self.pos_data)):
      #    print(f"{self.pos_data[i]}, {self.dose_data_raw[i]}, {self.dose_data[i]}")



   def output_plot(self):
      # build axis for plot:
      # x: positions
      # y: data, derivative, mean, inflection points
      y = self.derivative
      y2 = self.second_derivative
      y3 = [np.mean(self.dose_data) / 3.0 for _ in range(len(self.dose_data))]

      # plot components
      fig, ax = plt.subplots(2, 1)
      ax[0].plot(self.pos_data, self.dose_data, c = "blue")
      ax[0].scatter(
         [self.inflection_points[0][0], self.inflection_points[1][0], self.inflection_points[2][0], self.inflection_points[3][0], self.inflection_points[4][0], self.inflection_points[5][0]],
         [self.inflection_points[0][1][1], self.inflection_points[1][1][1], self.inflection_points[2][1][1], self.inflection_points[3][1][1], self.inflection_points[4][1][1], self.inflection_points[5][1][1]],
         c = "black")
      ax[1].plot(self.pos_data, self.derivative, c = "red")
      ax[1].plot(self.pos_data, self.second_derivative, c = "green")
      ax[1].scatter(
         [self.inflection_points[0][0], self.inflection_points[1][0], self.inflection_points[2][0], self.inflection_points[3][0], self.inflection_points[4][0], self.inflection_points[5][0]],
         [self.inflection_points[0][1][0], self.inflection_points[1][1][0], self.inflection_points[2][1][0], self.inflection_points[3][1][0], self.inflection_points[4][1][0], self.inflection_points[5][1][0]],
         c = "black")
      # save plot in the right profile subdirectory
      plt.savefig(f"{self.profile_out_dir}/{self.fields['SCAN_DEPTH']}.png")
      plt.close(fig)
      

   def produce_results(self):
      self.find_inflection_points()
      self.output_plot()




# --- Profile ---
# reads a file and extracts scan data present in that file
# by parsing scan by scan each field and collecting numerical data,
# organised in DataScans instances.

class Profile:
   name: str
   scans: list
   iso_field_size: str

   def __init__(self, file_path: str, out_dir: str, _test_settings: TestSettings):
      self.name = file_path.split("/")[-1].replace(".mcc", "")
      self.iso_field_size = self.name.split(" ")[3]
      self.scans = []
      
      raw_data = []

      print(f"{self.name}")

      # create output profile subdirectory
      # in the output directory there will be a subdir for every profile (with
      # the plot of each scan) and the inflection points table.
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

            self.scans.append(Scan(scan_num, raw_data, i, end_region_index, profile_out_dir, _test_settings))

            # skip content in between
            i = end_region_index
         i += 1



# --- Cacher ---
# stores all Profiles present in a folder

class Cacher:
   profiles: list
   res_dir: str # where to find data
   out_dir: str # where to store data

   # [!] filter faulty scans and log them to a file
   faulty_scans_list: list

   def __init__(self, _res_dir: str, _out_dir: str, _test_settings):
      self.profiles = []
      self.res_dir = _res_dir
      self.out_dir = _out_dir
      self.faulty_scans_list = []

      if not os.path.exists(_out_dir):
         os.mkdir(_out_dir)

      # create filelist
      filelist = []
      for root, dirs, files in os.walk(_res_dir):
         for _file in files:
            current_file_path = os.path.join(root,_file).replace("\\", "/")
            filelist.append(current_file_path)

      # create profile out of each file in filelist
      for f in filelist:
         self.profiles.append(Profile(f, self.out_dir, _test_settings))

      self.output_tables()
      


   # [!] transpose table to work with singular profiles rather than scans across profiles
   # [!] order profiles' scans and check for missing ones
   def create_table(self, profiles: list, measurement_depths: list):
      table = [["depth"] + [str(m / 10.0) for m in measurement_depths]]
      for p in range(len(profiles)):
         temp_table_row = [profiles[p].name.split(" ")[3]]
         for s in range(len(measurement_depths)):
            temp_profile_scans = profiles[p].scans
            dose_at_zero_index = temp_profile_scans[s].pos_data.index(0.0)
            dose_at_zero = temp_profile_scans[s].dose_data[dose_at_zero_index]
            if float(temp_profile_scans[s].fields["SCAN_DEPTH"]) == measurement_depths[s]:
               temp_table_row.append([
                  round(temp_profile_scans[s].inflection_points[0][0] / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[1][0] / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[0][1][1], 3),
                  round(temp_profile_scans[s].inflection_points[1][1][1], 3),
                  round(dose_at_zero, 3),
                  round(temp_profile_scans[s].inflection_points[2][0] / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[3][0] / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[2][1][1], 3),
                  round(temp_profile_scans[s].inflection_points[3][1][1], 3),
                  round(temp_profile_scans[s].inflection_points[4][0] / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[5][0] / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[4][1][1], 3),
                  round(temp_profile_scans[s].inflection_points[5][1][1], 3),
               ])
            else:
               temp_table_row.append("/")
               # insert null profile to offset right profile
               temp_profile_scans.insert(s, None)
         table.append(temp_table_row)
      return table

   def output_tables(self):
      # group profiles of same measurement
      measurements = dict()
      for p in self.profiles:
         measurement_name = p.name.split(" ")[-1]
         if measurement_name in measurements:
            measurements[measurement_name].append(p)
         else:
            measurements[measurement_name] = [p]

      for k, v in measurements.items():
         # loop through profiles to get all possible depth and order them
         measurement_depths = []
         for p in v:
            for s in p.scans:
               if not float(s.fields["SCAN_DEPTH"]) in measurement_depths:
                  measurement_depths.append(float(s.fields["SCAN_DEPTH"]))
         measurement_depths.sort()
         table = self.create_table(v, measurement_depths)
         table = transpose_table(table)
         with open(f"{self.out_dir}/{v[0].name.split(' ')[-1]}.txt", "w") as f:
            f.write("deriv_1, deriv_2, dose_1, dose_2, dose_0, second_deriv_1, second_deriv_2, second_dose_1, second_dose_2\n\n")
            for r in table:
               for c in r:
                  if isinstance(c, list):
                     str_list = f"{', '.join(map(str,c))}\t"
                     f.write(str_list.expandtabs(110))
                  else:
                     f.write(f"{c}\t".expandtabs(110))
               f.write("\n")
