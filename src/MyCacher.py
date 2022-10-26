# TODO: remove \t from Profile->raw_data, and avoid those horrible "\t\tSOMETHING"...

import os
import threading
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import MyUtils as utils

mpl.use('Agg')

# --- Scan ---
# stores scan data and parameters in a dictionary.


class Scan:
    # scan properties
    scan_num: int
    fields: dict
    pos_data: list
    dose_data: list
    init_bin: float
    fin_bin: float
    wanted_bin: float

    # calculated data
    first_derivative: list
    second_derivative: list
    third_derivative: list
    rebinned_pos_data: list
    rebinned_dose_data: list
    inflection_points: list
    lt25mm: bool
    dose_at_zero: float
    d1_left_fit_args: list
    d1_right_fit_args: list

    # useful data
    profile_out_dir: str

    def __init__(self, _scan_num, _raw_data, _begin, _end, _profile_out_dir, binning: float):
        self.scan_num = _scan_num
        self.fields = {}
        self.inflection_points = []
        self.lt25mm = False
        self.d1_left_fit_args: list
        self.d1_right_fit_args: list
        self.profile_out_dir = _profile_out_dir
        self.init_bin = 0.0
        self.fin_bin = 0.0
        self.wanted_bin = binning

        pos_data = []
        dose_data = []

        # ..:: parsing ::..
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
                    pos_data.append(float(nums[0]))
                    dose_data.append(float(nums[1]))
                    j += 1
                # stop parsing
                break
            i += 1

        dose_at_zero_index = pos_data.index(0.0)
        self.dose_at_zero = dose_data[dose_at_zero_index]

        # rebinning
        rebinned_pos_data = []
        rebinned_dose_data = []
        self.init_bin = abs(pos_data[0] - pos_data[1])
        # insert between every value and the next (if it exists) init_bin/wanted_bin number of interpolated points
        for i in range(len(pos_data)):
            inext = i + 1
            if inext == len(pos_data):
                break  # if there is no next point, stop before adding an excessive point
            # start interpolation
            t = 0.0
            while t < 1.0:
                rebinned_pos_data.append(utils.lerp(pos_data[i], pos_data[inext], t))
                rebinned_dose_data.append(
                    utils.lerp(dose_data[i], dose_data[inext], t))
                # increment by right amount
                t += (self.wanted_bin / self.init_bin)
        self.fin_bin = abs(rebinned_pos_data[0] - rebinned_pos_data[1])

        # ..:: fit ::..
        # peak position is gaussian center, and peak value is f(center)
        # with f beign gauss function with proper arguments for left and right fits
        first_derivative = utils.calc_derivative(pos_data, dose_data)
        self.orig_first_derivative = first_derivative[:]

        try:
            peak_pos = pos_data[first_derivative.index(max(first_derivative))]
            initial_parameters = [
                # [max(first_derivative), peak_pos, 10, 0],  # left fit
                # [min(first_derivative), -peak_pos, 10, 0],  # right fit
                [0.05, peak_pos, 3.0, -peak_pos / 2.0, max(first_derivative)],  # left fit
                [0.05, -peak_pos, 3.0, peak_pos / 2.0, max(first_derivative)],  # right fit
            ]
            # self.d1_left_fit_args, left_pcov = curve_fit(utils.gauss, pos_data[:len(
            #     pos_data) // 2], first_derivative[:len(first_derivative) // 2], initial_parameters[0])
            # self.d1_right_fit_args, right_pcov = curve_fit(utils.gauss, pos_data[len(
            #     pos_data) // 2:], first_derivative[len(first_derivative) // 2:], initial_parameters[1])
            self.d1_left_fit_args, left_pcov = curve_fit(utils.skew_normal_fit, pos_data[:len(pos_data) // 2], first_derivative[:len(first_derivative) // 2], initial_parameters[0])
            self.d1_right_fit_args, right_pcov = curve_fit(utils.skew_normal_fit, pos_data[len(pos_data) // 2:], first_derivative[len(first_derivative) // 2:], initial_parameters[1])
        except Exception as e:
            print(f"[error]: {e}")
            self.d1_left_fit_args = np.array([1, 1, 1, 1, 1])
            self.d1_right_fit_args = np.array([1, 1, 1, 1, 1])

        # half-bin correction
        self.d1_left_fit_args[1] += 0.25
        self.d1_right_fit_args[1] += 0.25

        # print(f"""{
        #          [list(map(round,initial_parameters[0],[5,5,5])), list(map(round,initial_parameters[1],[5,5,5]))]
        #       }\t{
        #          [list(map(round,self.d1_left_fit_args.tolist(),[5,5,5])), list(map(round,self.d1_right_fit_args.tolist(),[5,5,5]))]
        #       }""".expandtabs(8))

        # discretizing derivatives
        first_derivative = [utils.skew_normal(x, *self.d1_left_fit_args) for x in rebinned_pos_data[:len(rebinned_pos_data) // 2]]
        first_derivative.extend([utils.skew_normal(x, *self.d1_right_fit_args) for x in rebinned_pos_data[len(rebinned_pos_data) // 2:]])
        second_derivative = [utils.gauss_first_derivative( x, *self.d1_left_fit_args) for x in rebinned_pos_data[:len(rebinned_pos_data) // 2]]
        second_derivative.extend([utils.gauss_first_derivative( x, *self.d1_right_fit_args) for x in rebinned_pos_data[len(rebinned_pos_data) // 2:]])
        third_derivative = [utils.gauss_second_derivative( x, *self.d1_left_fit_args) for x in rebinned_pos_data[:len(rebinned_pos_data) // 2]]
        third_derivative.extend([utils.gauss_second_derivative( x, *self.d1_right_fit_args) for x in rebinned_pos_data[len(rebinned_pos_data) // 2:]])

        # ..:: inflection points ::..

        d1max1, d1maxi1, d1min1, d1mini1 = utils.max_and_min_in_range(
            first_derivative, None, None)

        d2max1, d2maxi1, d2min1, d2mini1 = utils.max_and_min_in_range(
            second_derivative, None, len(second_derivative) // 2)
        d2max2, d2maxi2, d2min2, d2mini2 = utils.max_and_min_in_range(
            second_derivative, len(second_derivative) // 2, None)

        d3max1, d3maxi1, d3min1, d3mini1 = utils.max_and_min_in_range(
            third_derivative, None, d2mini1)
        d3max2, d3maxi2, d3min2, d3mini2 = utils.max_and_min_in_range(
            third_derivative, d2mini1, d2mini2)
        d3max3, d3maxi3, d3min3, d3mini3 = utils.max_and_min_in_range(
            third_derivative, d2mini2, None)

        # additional point: dose(pos(d1max) - 25) exists ? eq25mm : lt25mm
        dose_offset_point_data = [0, [0, 0]]
        position_offset = 25.0
        try:
            dose_offset_point_data = [rebinned_pos_data[d1maxi1] - position_offset, [
                0, rebinned_dose_data[rebinned_pos_data.index(rebinned_pos_data[d1maxi1] - position_offset)]]]
        except:
            dose_offset_point_data = [
                rebinned_pos_data[0], [0, rebinned_dose_data[0]]]
            self.lt25mm = True

        # save inflection points data: [pos, [derivative_value, data_value]]
        self.inflection_points = [
            # d1
            [rebinned_pos_data[d1maxi1], [d1max1, rebinned_dose_data[d1maxi1]]],
            [rebinned_pos_data[d1mini1], [d1min1, rebinned_dose_data[d1mini1]]],
            # d2
            [rebinned_pos_data[d2maxi1], [d2max1, rebinned_dose_data[d2maxi1]]],
            [rebinned_pos_data[d2mini1], [d2min1, rebinned_dose_data[d2mini1]]],
            [rebinned_pos_data[d2mini2], [d2min2, rebinned_dose_data[d2mini2]]],
            [rebinned_pos_data[d2maxi2], [d2max2, rebinned_dose_data[d2maxi2]]],
            # d3
            [rebinned_pos_data[d3maxi1], [d3max1, rebinned_dose_data[d3maxi1]]],
            [rebinned_pos_data[d3mini1], [d3min1, rebinned_dose_data[d3mini1]]],
            [rebinned_pos_data[d3maxi2], [d3max2, rebinned_dose_data[d3maxi2]]],
            [rebinned_pos_data[d3mini2], [d3min2, rebinned_dose_data[d3mini2]]],
            [rebinned_pos_data[d3maxi3], [d3max3, rebinned_dose_data[d3maxi3]]],
            [rebinned_pos_data[d3mini3], [d3min3, rebinned_dose_data[d3mini3]]],
            # additional points
            dose_offset_point_data,
        ]

        # ..:: plotting ::..

        # params
        marker_size = 5
        line_width = 0.8

        # grouping inflection points for first and second scatter plots
        scatter_points_x = []
        scatter_points_y1 = []
        scatter_points_y2 = []
        for i in range(len(self.inflection_points)):
            scatter_points_x.append(self.inflection_points[i][0])
            scatter_points_y1.append(self.inflection_points[i][1][1])
            scatter_points_y2.append(self.inflection_points[i][1][0])

        fig, ax = plt.subplots(2, 1)

        ax[0].plot(pos_data, dose_data, c="blue", linewidth=line_width)
        ax[0].scatter(scatter_points_x, scatter_points_y1,
                      c="black", s=marker_size, zorder=9)
        ax[0].scatter(self.inflection_points[0][0], self.inflection_points[0]
                      [1][1], marker="+", c="cyan", zorder=10)
        ax[0].scatter(self.inflection_points[1][0], self.inflection_points[1]
                      [1][1], marker="+", c="cyan", zorder=10)

        ax[1].plot(pos_data, self.orig_first_derivative,
                   c="blue", linewidth=line_width)
        ax[1].plot(rebinned_pos_data, first_derivative,
                   c="red", linewidth=line_width)
        ax[1].plot(rebinned_pos_data, second_derivative,
                   c="green", linewidth=line_width)
        ax[1].plot(rebinned_pos_data, third_derivative,
                   c="purple", linewidth=line_width)
        ax[1].scatter(scatter_points_x, scatter_points_y2,
                      c="black", s=marker_size, zorder=9)
        ax[1].scatter(self.inflection_points[0][0], self.inflection_points[0]
                      [1][0], marker="+", c="cyan", zorder=10)
        ax[1].scatter(self.inflection_points[1][0], self.inflection_points[1]
                      [1][0], marker="+", c="cyan", zorder=10)

        # save plot in the right profile subdirectory
        plt.savefig(f"{self.profile_out_dir}/{self.fields['SCAN_DEPTH']}.png")
        plt.close(fig)

    def apply_filter(self, data, filter_name, arg):
        if filter_name == "moving_average":
            return utils.moving_average(data, arg)
        elif filter_name == "median_filter":
            return utils.median_filter(data, arg)
        # elif filter_name == "spline":
        #    tck = splrep(pos_data[:-1], dataset_index,)
        #    return splev(pos_data[:-1], tck).tolist()
        else:
            raise Exception("not yet implemented")


# --- Profile ---
# reads a file and extracts scan data present in that file
# by parsing scan by scan each field and collecting numerical data,
# organised in DataScans instances.

class Profile:
    name: str
    scans: list
    iso_field_size: str

    def __init__(self, file_path: str, out_dir: str, binning: float, shared_profiles: list):
        self.name = file_path.split("/")[-1].replace(".mcc", "")
        self.iso_field_size = self.name.split(" ")[3]
        self.scans = []

        raw_data = []

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
                    end_region_index = raw_data.index(
                        f"\tEND_SCAN  {scan_num}")
                else:
                    end_region_index = raw_data.index(f"\tEND_SCAN {scan_num}")

                self.scans.append(Scan(scan_num, raw_data, i,
                                  end_region_index, profile_out_dir, binning))

                # skip content in between
                i = end_region_index
            i += 1

        shared_profiles.append(self)
        print(f"finishing: {self.name}")


# --- Cacher ---
# stores all Profiles present in a folder

class Cacher:
    profiles: list
    res_dir: str  # where to find data
    out_dir: str  # where to store data

    def __init__(self, _res_dir: str, _out_dir: str, binning: float):
        self.profiles = []
        self.res_dir = _res_dir
        self.out_dir = _out_dir

        if not os.path.exists(_out_dir):
            os.mkdir(_out_dir)

        # create filelist
        filelist = []
        for root, dirs, files in os.walk(_res_dir):
            for _file in files:
                current_file_path = os.path.join(
                    root, _file).replace("\\", "/")
                filelist.append(current_file_path)
            break

        # create profile out of each file in filelist
        threads = []
        for f in filelist:
            t = threading.Thread(target=Profile, args=(
                f, self.out_dir, binning, self.profiles, ))
            print(f"starting: {f.split('/')[-1][:-4]}")
            t.start()
            threads.append(t)
            # self.profiles.append(Profile(f, self.out_dir, binning))
        for t in threads:
            t.join()

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
            field_size = int(
                round(float(profiles[p].scans[0].fields["FIELD_CROSSPLANE"]) / 100.0, 1) * 10.0)
            temp_table_row = [f"{field_size}X{field_size}"]
            # for each profile's scan
            for s in range(len(measurement_depths)):
                # some usefult values to remove noise in the code
                temp_profile_scans = profiles[p].scans

                # if this is the right measurement
                if float(temp_profile_scans[s].fields["SCAN_DEPTH"]) == measurement_depths[s]:

                    # fill the cell with dose_at_zero + inflection points data
                    # column_content = [round(dose_at_zero, 3)]
                    # for i in range(0, len(temp_profile_scans[s].inflection_points) - 1, 2):
                    #   # d1 positions and doses
                    #   column_content.append(round(temp_profile_scans[s].inflection_points[i][0] / 10.0, 3))
                    #   column_content.append(round(temp_profile_scans[s].inflection_points[i+1][0] / 10.0, 3))
                    #   column_content.append(round(temp_profile_scans[s].inflection_points[i][1][1], 3))
                    #   column_content.append(round(temp_profile_scans[s].inflection_points[i+1][1][1], 3))

                    temp_table_row.append([
                        # binning
                        round(temp_profile_scans[s].init_bin, data_precision),
                        round(temp_profile_scans[s].fin_bin, data_precision),
                        # D(0)
                        round(temp_profile_scans[s].dose_at_zero, data_precision),
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
                    # temp_table_row.append(column_content)

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
                        measurement_depths.append(
                            float(s.fields["SCAN_DEPTH"]))
            measurement_depths.sort()

            table = self.create_table(v, measurement_depths)

            # writing produced table to output file
            tabextension = 12
            with open(f"{self.out_dir}/{v[0].name.split(' ')[-1]}.txt", "w") as f:
                f.write("FS\td_cm\tini_bin\tbin\tD(0)\tp1d1sx\tD(p1d1sx)\tp1d1dx\tD(p1d1dx)\tp1d2sx\tp2d2sx\tD(p1d2sx)\tD(p2d2sx)\tp1d2dx\tp2d2dx\tD(p1d2dx)\tD(p2d2dx)\tp1d3sx\tp2d3sx\tp3d3sx\tD(p1d3sx)\tD(p2d3sx)\tD(p3d3sx)\tp1d3dx\tp2d3dx\tp3d3dx\tD(p1d3dx)\tD(p2d3dx)\tD(p3d3dx)\tD(IP-25)\toff25\tparams\n".expandtabs(tabextension))
                for row in range(1, len(table)):
                    for cell in range(1, len(table[row])):
                        cell_data = '\t'.join([str(x)
                                              for x in table[row][cell]])
                        f.write(
                            f"{table[row][0]}\t{table[0][cell]}\t{cell_data}\n".expandtabs(tabextension))

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
