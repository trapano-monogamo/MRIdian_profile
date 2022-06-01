# TODO: remove \t from Profile->raw_data, and avoid those horrible "\t\tSOMETHING"...

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from MyUtils import *




class ProcessingSettings:
   datasets: list
   filters:list
   args:list

   def __init__(self, _datasets:list, _filters:str, _args):
      self.datasets = _datasets
      self.filters = _filters
      self.args = _args


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
   lt25mm: bool

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
      self.lt25mm = False
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
            return moving_average(self.datasets[dataset_index], self.processing_settings.args[index_of_filter_name])
         elif self.processing_settings.filters[index_of_filter_name] == "median_filter":
            return median_filter(dataset_index, self.processing_settings.args[index_of_filter_name])
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
      # first derivative d1
      d1max = max(self.first_derivative)
      d1maxi = self.first_derivative.index(d1max)
      d1min = min(self.first_derivative)
      d1mini = self.first_derivative.index(d1min)

      # second derivative d2
      # 1
      d2max = max(self.second_derivative[: len(self.second_derivative) // 2])
      d2maxi = self.second_derivative.index(d2max, 0, len(self.second_derivative) // 2)
      d2min = min(self.second_derivative[: len(self.second_derivative) // 2])
      d2mini = self.second_derivative.index(d2min, 0, len(self.second_derivative) // 2)
      # 2
      d2max2 = max(self.second_derivative[len(self.second_derivative) // 2 :])
      d2maxi2 = self.second_derivative.index(d2max2, len(self.second_derivative) // 2)
      d2min2 = min(self.second_derivative[len(self.second_derivative) // 2 :])
      d2mini2 = self.second_derivative.index(d2min2, len(self.second_derivative) // 2)

      # third derivative d3
      # max
      d3max = max(self.third_derivative[:d1maxi])
      d3maxi = self.third_derivative.index(d3max, 0, d1maxi)
      d3max2 = max(self.third_derivative[d1maxi : len(self.third_derivative) // 2])
      d3maxi2 = self.third_derivative.index(d3max2, d1maxi, len(self.third_derivative) // 2)
      d3max3 = max(self.third_derivative[len(self.third_derivative) // 2 :])
      d3maxi3 = self.third_derivative.index(d3max3, len(self.third_derivative) // 2)
      # min
      d3min = min(self.third_derivative[:len(self.third_derivative) // 2])
      d3mini = self.third_derivative.index(d3min, 0, len(self.third_derivative) // 2)
      d3min2 = min(self.third_derivative[len(self.third_derivative) // 2 : d1mini])
      d3mini2 = self.third_derivative.index(d3min2, len(self.third_derivative) // 2, d1mini)
      d3min3 = min(self.third_derivative[d1mini:])
      d3mini3 = self.third_derivative.index(d3min3, d1mini)

      # dose(pos(d1max) - 30)
      dose_offset_point_data = [0, [0,0]]
      position_offset = 25.0
      try:
         dose_offset_point_data = [ self.pos_data[d1maxi] - position_offset, [0, self.dose_data[self.pos_data.index(self.pos_data[d1maxi] - position_offset)]] ]
      except:
         dose_offset_point_data = [ self.pos_data[0], [0, self.dose_data[0]] ]
         self.lt25mm = True

      # save inflection points data: [pos, [derivative_value, data_value]]
      self.inflection_points = [
         # d1
         [ (self.pos_data[d1maxi] + self.pos_data[d1maxi + 1]) / 2.0, [d1max, (self.dose_data[d1maxi] + self.dose_data[d1maxi+1]) / 2.0] ],
         [ (self.pos_data[d1mini] + self.pos_data[d1mini + 1]) / 2.0, [d1min, (self.dose_data[d1mini] + self.dose_data[d1mini+1]) / 2.0] ],
         # d2
         [ self.pos_data[d2maxi], [d2max, self.dose_data[d2maxi]] ],
         [ self.pos_data[d2mini], [d2min, self.dose_data[d2mini]] ],
         [ self.pos_data[d2mini2], [d2min2, self.dose_data[d2mini2]] ],
         [ self.pos_data[d2maxi2], [d2max2, self.dose_data[d2maxi2]] ],
         # d3
         [ self.pos_data[d3maxi], [d3max, self.dose_data[d3maxi]] ],
         [ self.pos_data[d3mini], [d3min, self.dose_data[d3mini]] ],
         [ self.pos_data[d3maxi2], [d3max2, self.dose_data[d3maxi2]] ],
         [ self.pos_data[d3mini2], [d3min2, self.dose_data[d3mini2]] ],
         [ self.pos_data[d3maxi3], [d3max3, self.dose_data[d3maxi3]] ],
         [ self.pos_data[d3mini3], [d3min3, self.dose_data[d3mini3]] ],
         # additional points
         dose_offset_point_data,
      ]

      # print(f"pos: {self.inflection_points[0][0] == -self.inflection_points[1][0]}, deriv: {self.inflection_points[0][1][0] == -self.inflection_points[1][1][0]}, dose: {self.inflection_points[0][1][1] == self.inflection_points[1][1][1]}, filt: {pmaxi > intersections[0]},{pmini < intersections[1]}")

      # for i in range(len(self.pos_data)):
      #    print(f"{self.pos_data[i]}, {self.dose_data_raw[i]}, {self.dose_data[i]}")



   def output_plot(self):
      marker_size = 7

      scatter_points_x = []
      scatter_points_y1 = []
      scatter_points_y2 = []
      for i in range(len(self.inflection_points)):
         scatter_points_x.append(self.inflection_points[i][0])
         scatter_points_y1.append(self.inflection_points[i][1][1])
         scatter_points_y2.append(self.inflection_points[i][1][0])

      # plot components
      fig, ax = plt.subplots(2, 1)

      ax[0].plot(self.pos_data, self.dose_data, c = "blue")
      ax[0].scatter(scatter_points_x, scatter_points_y1, c = "black", s = marker_size, zorder = 9)
      ax[0].scatter(self.inflection_points[0][0], self.inflection_points[0][1][1], marker = "+", c = "cyan", zorder = 10)
      ax[0].scatter(self.inflection_points[1][0], self.inflection_points[1][1][1], marker = "+", c = "cyan", zorder = 10)

      ax[1].plot(self.pos_data, self.first_derivative, c = "red")
      ax[1].plot(self.pos_data, self.second_derivative, c = "green")
      ax[1].plot(self.pos_data, self.third_derivative, c = "purple")
      ax[1].scatter(scatter_points_x, scatter_points_y2, c = "black", s = marker_size, zorder = 9)
      ax[1].scatter(self.inflection_points[0][0], self.inflection_points[0][1][0], marker = "+", c = "cyan", zorder = 10)
      ax[1].scatter(self.inflection_points[1][0], self.inflection_points[1][1][0], marker = "+", c = "cyan", zorder = 10)

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
      


   def create_table(self, profiles: list, measurement_depths: list):
      # create table's first (transposed first column)
      table = [["depth"] + [str(m / 10.0) for m in measurement_depths]]
      # for each profile
      for p in range(len(profiles)):
         # prepare the table's row
         temp_table_row = [profiles[p].name.split(" ")[3]]
         # for each profile's scan
         for s in range(len(measurement_depths)):
            # some usefult values to remove noise in the code
            temp_profile_scans = profiles[p].scans
            dose_at_zero_index = temp_profile_scans[s].pos_data.index(0.0)
            dose_at_zero = temp_profile_scans[s].dose_data[dose_at_zero_index]

            # if this is the right measurement
            if float(temp_profile_scans[s].fields["SCAN_DEPTH"]) == measurement_depths[s]:

               # fill the cell with dose_at_zero + inflection points data
               column_content = [round(dose_at_zero, 3)]
               #for i in range(0, len(temp_profile_scans[s].inflection_points) - 1, 2):
               #   # d1 positions and doses
               #   column_content.append(round(temp_profile_scans[s].inflection_points[i][0] / 10.0, 3))
               #   column_content.append(round(temp_profile_scans[s].inflection_points[i+1][0] / 10.0, 3))
               #   column_content.append(round(temp_profile_scans[s].inflection_points[i][1][1], 3))
               #   column_content.append(round(temp_profile_scans[s].inflection_points[i+1][1][1], 3))

               temp_table_row.append([
                  round(dose_at_zero, 3),
                  # d1: max
                  round(temp_profile_scans[s].inflection_points[0][0]     / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[0][1][1]  / 10.0, 3),
                  # d1: min
                  round(temp_profile_scans[s].inflection_points[1][0]     / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[1][1][1]  / 10.0, 3),
                  # d2: max1, min1
                  round(temp_profile_scans[s].inflection_points[2][0]     / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[3][0]     / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[2][1][1]  / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[3][1][1]  / 10.0, 3),
                  # d2: min2, max2
                  round(temp_profile_scans[s].inflection_points[4][0]     / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[5][0]     / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[4][1][1]  / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[5][1][1]  / 10.0, 3),
                  # d3: max1, min1, max2
                  round(temp_profile_scans[s].inflection_points[6][0]     / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[7][0]     / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[8][0]     / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[6][1][1]  / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[7][1][1]  / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[8][1][1]  / 10.0, 3),
                  # d3: min2, max3, min3
                  round(temp_profile_scans[s].inflection_points[9][0]     / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[10][0]    / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[11][0]    / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[9][1][1]  / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[10][1][1] / 10.0, 3),
                  round(temp_profile_scans[s].inflection_points[11][1][1] / 10.0, 3),
                  # additional points
                  round(temp_profile_scans[s].inflection_points[-1][1][1] / 10.0, 3),
                  "lt25" if temp_profile_scans[s].lt25mm else "eq25",
               ])

               # append the cell to the row
               #temp_table_row.append(column_content)

            # if it's not the right scan
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

         # table stuff
         table = self.create_table(v, measurement_depths)

         tabextension = 12
         with open(f"{self.out_dir}/{v[0].name.split(' ')[-1]}.txt", "w") as f:
            f.write("FS\td_cm\tD(0)\tp1d1sx\tD(p1d1sx)\tp1d1dx\tD(p1d1dx)\tp1d2sx\tp2d2sx\tD(p1d2sx)\tD(p2d2sx)\tp1d2dx\tp2d2dx\tD(p1d2dx)\tD(p2d2dx)\tp1d3sx\tp2d3sx\tp3d3sx\tD(p1d3sx)\tD(p2d3sx)\tD(p3d3sx)\tp1d3dx\tp2d3dx\tp3d3dx\tD(p1d3dx)\tD(p2d3dx)\tD(p3d3dx)\tD(IP-25)\toff25\n".expandtabs(tabextension))
            for row in range(1, len(table)):
               for cell in range(1, len(table[row])):
                  cell_data = '\t'.join([str(x) for x in table[row][cell]])
                  f.write(f"{table[row][0]}\t{table[0][cell]}\t{cell_data}\n".expandtabs(tabextension))

         # table = transpose_table(table)
         # with open(f"{self.out_dir}/{v[0].name.split(' ')[-1]}.txt", "w") as f:
         #    tab_expansion = 450
         #    f.write("pos_d1_1, pos_d1_2, dose_d1_1, dose_d2_2, pos_d2_1, pos_d2_2, pos_d3_1, pos_d3_2, pos_d3_3, pos_d3_4, dose_0\n\n")
         #    for r in table:
         #       for c in r:
         #          if isinstance(c, list):
         #             str_list = f"{', '.join(map(str,c))}\t"
         #             f.write(str_list.expandtabs(tab_expansion))
         #          else:
         #             f.write(f"{c}\t".expandtabs(tab_expansion))
         #       f.write("\n")
