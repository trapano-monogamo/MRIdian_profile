from MyCacher import *
import multiprocessing
import time

def main():
   start_time = time.time()

   p1 = multiprocessing.Process(target = Cacher, args = ("./res/MR_LINAC/", f"./out/MR_LINAC/multithreaded test 0.1mm/", 0.1))
   p1.start()
   p2 = multiprocessing.Process(target = Cacher, args = ("./res/TRUEBEAM/", "./out/TRUEBEAM/multithreaded test 0.1mm", 0.1))
   p2.start()

   p1.join()
   p2.join()

   end_time = time.time()
   print(f"exec time: {end_time - start_time}")

if __name__ == '__main__':
   main()
