# TODO: remove \t from Profile->raw_data, and avoid those horrible "\t\tSOMETHING"...

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from MyUtils import *




class ProcessingSettings:
   datasets: list
   filters:list

   def __init__(self, _datasets:list, _filters:str):
      self.datasets = _datasets
      self.filters = _filters


# --- Scan ---
# stores scan data and parameters in a dictionary.

class Scan:
   # scan properties
   scan_num: int
   fields: dict
   pos_data: list
   dose_data: list

   # calculated data
   first_derivative: list
   second_derivative: list
   third_derivative: list
   inflection_points: list

   # useful data
   profile_out_dir: str

   # test driven stuff babyyy
   processing_settings: ProcessingSettings
   datasets: list

   def __init__(self, _scan_num, _raw_data, _begin, _end, _profile_out_dir, _processing_settings: ProcessingSettings):
      self.scan_num = _scan_num
      self.fields = {}
      self.pos_data = []
      self.dose_data = []
      self.first_derivative = []
      self.second_derivative = []
      self.third_derivative = []
      self.inflection_points = []
      self.profile_out_dir = _profile_out_dir

      self.datasets = [
         self.dose_data,
         self.first_derivative,
         self.second_derivative,
         self.third_derivative
      ]
      self.processing_settings = _processing_settings

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


   def apply_filter(self, dataset_index):
      if dataset_index in self.processing_settings.datasets:
         index_of_filter_name = self.processing_settings.datasets.index(dataset_index)
         if self.processing_settings.filters[index_of_filter_name] == "moving_average":
            return moving_average(self.datasets[dataset_index], 3)
         elif self.processing_settings.filters[index_of_filter_name] == "median_filter":
            return median_filter(dataset_index, 3)
         elif self.processing_settings.filters[index_of_filter_name] == "spline":
            tck = splrep(self.pos_data[:-1], dataset_index,)
            return splev(self.pos_data[:-1], tck).tolist()
         else:
            raise Exception("not yet implemented")
      else:
         return self.datasets[dataset_index]


   # eww... and it's not the emacs' browser
   def process_data(self):
      # dose data
      self.dose_data = self.apply_filter(self.datasets.index(self.dose_data))
      self.datasets[0] = self.first_derivative

      # derivatives
      # first derivative and filter
      self.first_derivative = calc_derivative(self.pos_data, self.dose_data)
      self.datasets[1] = self.first_derivative
      self.first_derivative = self.apply_filter(self.datasets.index(self.first_derivative))
      # second derivative and filter (and correction)
      self.second_derivative = calc_derivative(self.pos_data, self.first_derivative)
      self.datasets[2] = self.second_derivative
      self.second_derivative = self.apply_filter(self.datasets.index(self.second_derivative))
      ranges = [[0,10], [len(self.second_derivative) - 10, len(self.second_derivative)]]
      self.second_derivative = median_filter(self.second_derivative, 3, ranges)
      # third derivative
      self.third_derivative = calc_derivative(self.pos_data, self.second_derivative)
      self.datasets[3] = self.third_derivative
      self.third_derivative = self.apply_filter(self.datasets.index(self.third_derivative))

      # other data
      # -


   # Calculate derivative, than apply median filter to portions of the derivative to
   # eliminate peaks generated by experimental errors jumps in the original data.
   # These portions are identified as external intervals of the intersections between
   # the original data curve and a third of the mean (arbitrary value... ugh).
   # Once filtered, max and min of the derivative will be the functions high and low,
   # which are the derivative values at the inflection points positions.
   # To find the positions take the index of the peak in the derivative, and use it as index for
   # both derivative function and original data
   def find_inflection_points(self):
      self.process_data()

      # find peaks/valleys and their positions
      pmax = max(self.first_derivative)
      pmaxi = self.first_derivative.index(pmax)
      pmin = min(self.first_derivative)
      pmini = self.first_derivative.index(pmin)

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
      y = self.first_derivative
      y2 = self.second_derivative
      y3 = [np.mean(self.dose_data) / 3.0 for _ in range(len(self.dose_data))]

      # plot components
      fig, ax = plt.subplots(2, 1)
      ax[0].plot(self.pos_data, self.dose_data, c = "blue")
      ax[0].scatter(
         [self.inflection_points[0][0], self.inflection_points[1][0], self.inflection_points[2][0], self.inflection_points[3][0], self.inflection_points[4][0], self.inflection_points[5][0]],
         [self.inflection_points[0][1][1], self.inflection_points[1][1][1], self.inflection_points[2][1][1], self.inflection_points[3][1][1], self.inflection_points[4][1][1], self.inflection_points[5][1][1]],
         c = "black")
      ax[1].plot(self.pos_data, self.first_derivative, c = "red")
      ax[1].plot(self.pos_data, self.second_derivative, c = "green")
      ax[1].plot(self.pos_data, self.third_derivative, c = "purple")
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

   def __init__(self, file_path: str, out_dir: str, _processing_settings: ProcessingSettings):
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

            self.scans.append(Scan(scan_num, raw_data, i, end_region_index, profile_out_dir, _processing_settings))

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

   def __init__(self, _res_dir: str, _out_dir: str, _processing_settings):
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
         self.profiles.append(Profile(f, self.out_dir, _processing_settings))

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
