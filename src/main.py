from MyCacher import *
import multiprocessing

def main():
   p1 = multiprocessing.Process(target = Cacher, args = ("./res/MR_LINAC/", f"./out/MR_LINAC/multithreaded test #2 0.001mm/", 0.001))
   p1.start()
   p2 = multiprocessing.Process(target = Cacher, args = ("./res/TRUEBEAM/", "./out/TRUEBEAM/multithreaded test #2 0.001mm", 0.001))
   p2.start()

   p1.join()
   p2.join()

if __name__ == '__main__':
   main()
