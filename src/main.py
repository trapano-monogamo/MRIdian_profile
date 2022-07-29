from MyCacher import *
import time

def main():
   start_time = time.time()

   Cacher(
      "./res/MR_LINAC/",
      f"./out/MR_LINAC/singlethreaded test 0.01mm/",
      0.01)
   Cacher(
      "./res/TRUEBEAM/",
      "./out/TRUEBEAM/singlethreaded test 0.01mm",
      0.01)

   end_time = time.time()
   print(f"exec time: {end_time - start_time}")

if __name__ == '__main__':
   main()
