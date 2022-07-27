# TODO: remove \t from Profile->raw_data, and avoid those horrible "\t\tSOMETHING"...

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lmfit.models import GaussianModel
from MyUtils import *




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
   d1_left_fit_args: list
   d1_right_fit_args: list
   binning: float

   # useful data
   profile_out_dir: str

   def __init__(self, _scan_num, _raw_data, _begin, _end, _profile_out_dir, binning:float):
      self.scan_num = _scan_num
      self.fields = {}
      self.pos_data = []
      self.dose_data = []
      self.first_derivative = []
      self.second_derivative = []
      self.third_derivative = []
      self.inflection_points = []
      self.lt25mm = False
      self.d1_left_fit_args: list
      self.d1_right_fit_args: list
      self.profile_out_dir = _profile_out_dir
      self.binning = binning

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


   def apply_filter(self, data, filter_name, arg):
      if filter_name == "moving_average":
         return moving_average(data, arg)
      elif filter_name == "median_filter":
         return median_filter(data, arg)
      # elif filter_name == "spline":
      #    tck = splrep(self.pos_data[:-1], dataset_index,)
      #    return splev(self.pos_data[:-1], tck).tolist()
      else:
         raise Exception("not yet implemented")


   def process_data(self):
      self.rebinned_pos_data = []
      self.rebinned_dose_data = []

      # rebinning
      current_binning = abs(self.pos_data[0] - self.pos_data[1])
      for i in range(len(self.pos_data)):
         inext = i + 1
         if inext == len(self.pos_data): break
         t = 0.0
         while t < 1.0:
            self.rebinned_pos_data.append(lerp(self.pos_data[i], self.pos_data[inext], t))
            self.rebinned_dose_data.append(lerp(self.dose_data[i], self.dose_data[inext], t))
            t += (self.binning / current_binning)

      self.first_derivative = calc_derivative(self.pos_data, self.dose_data)
      self.orig_first_derivative = self.first_derivative[:]
      try:
         initial_parameters = [
            [max(self.first_derivative), -(abs(self.pos_data[0]) // 5 * 4), 20], # left fit
            [min(self.first_derivative), self.pos_data[-1] // 5 * 4, 20], # right fit
         ]
         self.d1_left_fit_args, left_pcov = curve_fit(gauss, self.pos_data[:len(self.pos_data) // 2], self.first_derivative[:len(self.first_derivative) // 2], initial_parameters[0])
         self.d1_right_fit_args, right_pcov = curve_fit(gauss, self.pos_data[len(self.pos_data) // 2:], self.first_derivative[len(self.first_derivative) // 2:], initial_parameters[1])
      except:
         self.d1_left_fit_args = np.array([1,1,1])
         self.d1_right_fit_args = np.array([1,1,1])
      self.d1_left_fit_args[1] += 0.25
      self.d1_right_fit_args[1] += 0.25
      self.first_derivative = [gauss(x, *self.d1_left_fit_args) for x in self.rebinned_pos_data[:len(self.rebinned_pos_data) // 2]]
      self.first_derivative.extend([gauss(x, *self.d1_right_fit_args) for x in self.rebinned_pos_data[len(self.rebinned_pos_data) // 2:]])

      # !!!
      # peak position is gaussian center, and peak value is f(center)
      # with f beign gauss function with proper arguments for left and right fits

      # self.second_derivative = calc_derivative(self.pos_data, self.first_derivative)
      self.second_derivative = [gauss_first_derivative(x, *self.d1_left_fit_args) for x in self.rebinned_pos_data[:len(self.rebinned_pos_data)//2]]
      self.second_derivative.extend([gauss_first_derivative(x, *self.d1_right_fit_args) for x in self.rebinned_pos_data[len(self.rebinned_pos_data)//2:]])

      # self.third_derivative = calc_derivative(self.pos_data, self.second_derivative)
      self.third_derivative = [gauss_second_derivative(x, *self.d1_left_fit_args) for x in self.rebinned_pos_data[:len(self.rebinned_pos_data)//2]]
      self.third_derivative.extend([gauss_second_derivative(x, *self.d1_right_fit_args) for x in self.rebinned_pos_data[len(self.rebinned_pos_data)//2:]])

   def find_inflection_points(self):
      self.process_data()

      d1max1, d1maxi1, d1min1, d1mini1 = max_and_min_in_range(self.first_derivative, None, None)

      d2max1, d2maxi1, d2min1, d2mini1 = max_and_min_in_range(self.second_derivative, None, len(self.second_derivative) // 2)
      d2max2, d2maxi2, d2min2, d2mini2 = max_and_min_in_range(self.second_derivative, len(self.second_derivative) // 2, None)

      d3max1, d3maxi1, d3min1, d3mini1 = max_and_min_in_range(self.third_derivative, None, d2mini1)
      d3max2, d3maxi2, d3min2, d3mini2 = max_and_min_in_range(self.third_derivative, d2mini1, d2mini2)
      d3max3, d3maxi3, d3min3, d3mini3 = max_and_min_in_range(self.third_derivative, d2mini2, None)

      # dose(pos(d1max) - 30)
      dose_offset_point_data = [0, [0,0]]
      position_offset = 25.0
      try:
         dose_offset_point_data = [ self.rebinned_pos_data[d1maxi1] - position_offset, [0, self.rebinned_dose_data[self.rebinned_pos_data.index(self.rebinned_pos_data[d1maxi1] - position_offset)]] ]
      except:
         dose_offset_point_data = [ self.rebinned_pos_data[0], [0, self.rebinned_dose_data[0]] ]
         self.lt25mm = True


      # save inflection points data: [pos, [derivative_value, data_value]]
      # self.inflection_points = [
      #    # d1
      #    [ (self.rebinned_pos_data[d1maxi1] + self.rebinned_pos_data[d1maxi1 + 1]) / 2.0, [d1max1, (self.rebinned_dose_data[d1maxi1] + self.rebinned_dose_data[d1maxi1+1]) / 2.0] ],
      #    [ (self.rebinned_pos_data[d1mini1] + self.rebinned_pos_data[d1mini1 + 1]) / 2.0, [d1min1, (self.rebinned_dose_data[d1mini1] + self.rebinned_dose_data[d1mini1+1]) / 2.0] ],
      #    # d2
      #    [ (self.rebinned_pos_data[d2maxi1] + self.rebinned_pos_data[d2maxi1 + 1]) / 2.0, [d2max1, (self.rebinned_dose_data[d2maxi1] + self.rebinned_dose_data[d2maxi1 + 1]) / 2.0] ],
      #    [ (self.rebinned_pos_data[d2mini1] + self.rebinned_pos_data[d2mini1 + 1]) / 2.0, [d2min1, (self.rebinned_dose_data[d2mini1] + self.rebinned_dose_data[d2mini1 + 1]) / 2.0] ],
      #    [ (self.rebinned_pos_data[d2mini2] + self.rebinned_pos_data[d2mini2 + 1]) / 2.0, [d2min2, (self.rebinned_dose_data[d2mini2] + self.rebinned_dose_data[d2mini2 + 1]) / 2.0] ],
      #    [ (self.rebinned_pos_data[d2maxi2] + self.rebinned_pos_data[d2maxi2 + 1]) / 2.0, [d2max2, (self.rebinned_dose_data[d2maxi2] + self.rebinned_dose_data[d2maxi2 + 1]) / 2.0] ],
      #    # d3
      #    [ (self.rebinned_pos_data[d3maxi1] + self.rebinned_pos_data[d3maxi1 + 1]) / 2.0, [d3max1, (self.rebinned_dose_data[d3maxi1] + self.rebinned_dose_data[d3maxi1 + 1]) / 2.0] ],
      #    [ (self.rebinned_pos_data[d3mini1] + self.rebinned_pos_data[d3mini1 + 1]) / 2.0, [d3min1, (self.rebinned_dose_data[d3mini1] + self.rebinned_dose_data[d3mini1 + 1]) / 2.0] ],
      #    [ (self.rebinned_pos_data[d3maxi2] + self.rebinned_pos_data[d3maxi2 + 1]) / 2.0, [d3max2, (self.rebinned_dose_data[d3maxi2] + self.rebinned_dose_data[d3maxi2 + 1]) / 2.0] ],
      #    [ (self.rebinned_pos_data[d3mini2] + self.rebinned_pos_data[d3mini2 + 1]) / 2.0, [d3min2, (self.rebinned_dose_data[d3mini2] + self.rebinned_dose_data[d3mini2 + 1]) / 2.0] ],
      #    [ (self.rebinned_pos_data[d3maxi3] + self.rebinned_pos_data[d3maxi3 + 1]) / 2.0, [d3max3, (self.rebinned_dose_data[d3maxi3] + self.rebinned_dose_data[d3maxi3 + 1]) / 2.0] ],
      #    [ (self.rebinned_pos_data[d3mini3] + self.rebinned_pos_data[d3mini3 + 1]) / 2.0, [d3min3, (self.rebinned_dose_data[d3mini3] + self.rebinned_dose_data[d3mini3 + 1]) / 2.0] ],
      #    # additional points
      #    dose_offset_point_data,
      # ]
      self.inflection_points = [
         # d1
         [ self.rebinned_pos_data[d1maxi1], [d1max1, self.rebinned_dose_data[d1maxi1]] ],
         [ self.rebinned_pos_data[d1mini1], [d1min1, self.rebinned_dose_data[d1mini1]] ],
         # d2
         [ self.rebinned_pos_data[d2maxi1], [d2max1, self.rebinned_dose_data[d2maxi1]] ],
         [ self.rebinned_pos_data[d2mini1], [d2min1, self.rebinned_dose_data[d2mini1]] ],
         [ self.rebinned_pos_data[d2mini2], [d2min2, self.rebinned_dose_data[d2mini2]] ],
         [ self.rebinned_pos_data[d2maxi2], [d2max2, self.rebinned_dose_data[d2maxi2]] ],
         # d3
         [ self.rebinned_pos_data[d3maxi1], [d3max1, self.rebinned_dose_data[d3maxi1]] ],
         [ self.rebinned_pos_data[d3mini1], [d3min1, self.rebinned_dose_data[d3mini1]] ],
         [ self.rebinned_pos_data[d3maxi2], [d3max2, self.rebinned_dose_data[d3maxi2]] ],
         [ self.rebinned_pos_data[d3mini2], [d3min2, self.rebinned_dose_data[d3mini2]] ],
         [ self.rebinned_pos_data[d3maxi3], [d3max3, self.rebinned_dose_data[d3maxi3]] ],
         [ self.rebinned_pos_data[d3mini3], [d3min3, self.rebinned_dose_data[d3mini3]] ],
         # additional points
         dose_offset_point_data,
      ]

      # print(f"pos: {self.inflection_points[0][0] == -self.inflection_points[1][0]}, deriv: {self.inflection_points[0][1][0] == -self.inflection_points[1][1][0]}, dose: {self.inflection_points[0][1][1] == self.inflection_points[1][1][1]}, filt: {pmaxi > intersections[0]},{pmini < intersections[1]}")

      # for i in range(len(self.pos_data)):
      #    print(f"{self.pos_data[i]}, {self.dose_data_raw[i]}, {self.dose_data[i]}")



   def output_plot(self):
      marker_size = 5
      line_width = 0.8

      scatter_points_x = []
      scatter_points_y1 = []
      scatter_points_y2 = []
      for i in range(len(self.inflection_points)):
         scatter_points_x.append(self.inflection_points[i][0])
         scatter_points_y1.append(self.inflection_points[i][1][1])
         scatter_points_y2.append(self.inflection_points[i][1][0])

      # plot components
      fig, ax = plt.subplots(2, 1)

      ax[0].plot(self.pos_data, self.dose_data, c = "blue", linewidth = line_width)
      ax[0].scatter(scatter_points_x, scatter_points_y1, c = "black", s = marker_size, zorder = 9)
      ax[0].scatter(self.inflection_points[0][0], self.inflection_points[0][1][1], marker = "+", c = "cyan", zorder = 10)
      ax[0].scatter(self.inflection_points[1][0], self.inflection_points[1][1][1], marker = "+", c = "cyan", zorder = 10)

      ax[1].plot(self.pos_data, self.orig_first_derivative, c = "blue", linewidth = line_width)
      ax[1].plot(self.rebinned_pos_data, self.first_derivative, c = "red", linewidth = line_width)
      ax[1].plot(self.rebinned_pos_data, self.second_derivative, c = "green", linewidth = line_width)
      ax[1].plot(self.rebinned_pos_data, self.third_derivative, c = "purple", linewidth = line_width)
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

   def __init__(self, file_path: str, out_dir: str, binning: float):
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

            self.scans.append(Scan(scan_num, raw_data, i, end_region_index, profile_out_dir, binning))

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

   def __init__(self, _res_dir: str, _out_dir: str, binning: float):
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
         self.profiles.append(Profile(f, self.out_dir, binning))

      self.output_tables()
      


   def create_table(self, profiles: list, measurement_depths: list):
      # create table's first (transposed first column)
      table = [["depth"] + [str(m / 10.0) for m in measurement_depths]]
      data_precision = 5
      # for each profile
      for p in range(len(profiles)):
         # prepare the table's row
         # field size
         # temp_table_row = [profiles[p].name.split(" ")[3]]
         field_size = int(round(float(profiles[p].scans[0].fields["FIELD_CROSSPLANE"]) / 100.0, 1) * 10.0)
         temp_table_row = [f"{field_size}X{field_size}"]
         # for each profile's scan
         for s in range(len(measurement_depths)):
            # some usefult values to remove noise in the code
            temp_profile_scans = profiles[p].scans
            dose_at_zero_index = temp_profile_scans[s].pos_data.index(0.0)
            dose_at_zero = temp_profile_scans[s].dose_data[dose_at_zero_index]

            # if this is the right measurement
            if float(temp_profile_scans[s].fields["SCAN_DEPTH"]) == measurement_depths[s]:

               # fill the cell with dose_at_zero + inflection points data
               #column_content = [round(dose_at_zero, 3)]
               #for i in range(0, len(temp_profile_scans[s].inflection_points) - 1, 2):
               #   # d1 positions and doses
               #   column_content.append(round(temp_profile_scans[s].inflection_points[i][0] / 10.0, 3))
               #   column_content.append(round(temp_profile_scans[s].inflection_points[i+1][0] / 10.0, 3))
               #   column_content.append(round(temp_profile_scans[s].inflection_points[i][1][1], 3))
               #   column_content.append(round(temp_profile_scans[s].inflection_points[i+1][1][1], 3))

               temp_table_row.append([
                  # binning
                  round(abs(temp_profile_scans[s].pos_data[0] - temp_profile_scans[s].pos_data[1]), data_precision),
                  round(abs(temp_profile_scans[s].rebinned_pos_data[0] - temp_profile_scans[s].rebinned_pos_data[1]), data_precision),
                  # D(0)
                  round(dose_at_zero, data_precision),
                  # d1: max
                  round(temp_profile_scans[s].inflection_points[0][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[0][1][1], data_precision),
                  # d1: min
                  round(temp_profile_scans[s].inflection_points[1][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[1][1][1], data_precision),
                  # d2: max1, min1
                  round(temp_profile_scans[s].inflection_points[2][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[3][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[2][1][1], data_precision),
                  round(temp_profile_scans[s].inflection_points[3][1][1], data_precision),
                  # d2: min2, max2
                  round(temp_profile_scans[s].inflection_points[4][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[5][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[4][1][1], data_precision),
                  round(temp_profile_scans[s].inflection_points[5][1][1], data_precision),
                  # d3: max1, min1, max2
                  round(temp_profile_scans[s].inflection_points[6][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[7][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[8][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[6][1][1], data_precision),
                  round(temp_profile_scans[s].inflection_points[7][1][1], data_precision),
                  round(temp_profile_scans[s].inflection_points[8][1][1], data_precision),
                  # d3: min2, max3, min3
                  round(temp_profile_scans[s].inflection_points[9][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[10][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[11][0], data_precision),
                  round(temp_profile_scans[s].inflection_points[9][1][1], data_precision),
                  round(temp_profile_scans[s].inflection_points[10][1][1], data_precision),
                  round(temp_profile_scans[s].inflection_points[11][1][1], data_precision),
                  # additional points
                  round(temp_profile_scans[s].inflection_points[-1][1][1] / 10.0, data_precision),
                  "lt25" if temp_profile_scans[s].lt25mm else "eq25",
                  *[round(n, data_precision) for n in temp_profile_scans[s].d1_left_fit_args.tolist()
                     + temp_profile_scans[s].d1_right_fit_args.tolist()],
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

         table = self.create_table(v, measurement_depths)

         # writing produced table to output file
         tabextension = 12
         with open(f"{self.out_dir}/{v[0].name.split(' ')[-1]}.txt", "w") as f:
            f.write("FS\td_cm\tini_bin\tbin\tD(0)\tp1d1sx\tD(p1d1sx)\tp1d1dx\tD(p1d1dx)\tp1d2sx\tp2d2sx\tD(p1d2sx)\tD(p2d2sx)\tp1d2dx\tp2d2dx\tD(p1d2dx)\tD(p2d2dx)\tp1d3sx\tp2d3sx\tp3d3sx\tD(p1d3sx)\tD(p2d3sx)\tD(p3d3sx)\tp1d3dx\tp2d3dx\tp3d3dx\tD(p1d3dx)\tD(p2d3dx)\tD(p3d3dx)\tD(IP-25)\toff25\tparams\n".expandtabs(tabextension))
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
