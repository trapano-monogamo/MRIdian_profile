# --- DataScan ---
# stores scan data and parameters in a dictionary.

class DataScan:
   task_name: str
   scan_num: int
   data: dict
   begin: int
   end: int

   def __init__(self, _scan_num, _file, _begin, _end):
      self.scan_num = _scan_num
      self.begin = _begin
      self.end = _end

      # handmade parsing rules

   def log(self):
      print(f"n: {self.scan_num} , begin: {self.begin} , end: {self.end}")



# --- DataProfile ---
# reads a file and extracts scan data present in that file
# by parsing scan by scan each field and collecting numerical data,
# organised in DataScans instances.

class DataProfile:
   raw_data: list
   scans: list

   def __init__(self, file_path: str):
      self.scans = []
      self.raw_data = []

      # read file
      with open(file_path, "r") as f:
         self.raw_data = f.read().split("\n")
         
      # handmade parsing rules
      i = 0
      end_region_index = 0
      while i < len(self.raw_data):
         # when a BEGIN_SCAN is found, highlight its region up to END_SCAN, and build a DataScan out of it
         if "\tBEGIN_SCAN" in self.raw_data[i].split(' '):
            scan_num = int(self.raw_data[i].split(" ")[-1])

            # OMG!!! what are we, babies?! c'mon...
            if scan_num < 10:
               end_region_index = self.raw_data.index(f"\tEND_SCAN  {scan_num}")
            else:
               end_region_index = self.raw_data.index(f"\tEND_SCAN {scan_num}")

            self.scans.append(DataScan(scan_num, file_path, i, end_region_index))

            # skip content in between
            i = end_region_index
         i += 1

      for s in self.scans:
         s.log()

      # storing data



# --- DataCacher ---
# stores all DataProfiles present in a folder

class DataCacher:
   profiles: list

   def __init__(self):
      pass






def main():
   d = DataProfile("./res/X06 FFF OPEN 15X15 X Average_Symmetrised.mcc")
   d = DataProfile("./res/X06 FFF OPEN 15X15 X OLBIA_Symmetrised.mcc")

if __name__ == '__main__':
   main()
