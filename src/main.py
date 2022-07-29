from MyCacher import *
import multiprocessing
import time

def main():
   time_start = time.time()

   p1 = multiprocessing.Process(target = Cacher, args = ("./res/MR_LINAC/", f"./out/MR_LINAC/multithreaded test 0.01mm/", 0.01))
   p1.start()
   p2 = multiprocessing.Process(target = Cacher, args = ("./res/TRUEBEAM/", "./out/TRUEBEAM/multithreaded test 0.01mm", 0.01))
   p2.start()

   p1.join()
   p2.join()

   time_end = time.time()
   print(f"exec time: {time_end - time_start}")

if __name__ == '__main__':
   main()
