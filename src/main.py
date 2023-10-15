from MyCacher import *
from chi_adder import *
import multiprocessing
import time

def main():
    data_package = "MR_LINAC"

    # start_time = time.time()

    # p1 = multiprocessing.Process(target=Cacher, args=(f"./res/{data_package}/", f"./out/{data_package}/chi tabs/", 0.001))
    # p1.start()

    # p1.join()

    # end_time = time.time()
    # print(f"exec time: {end_time - start_time}")

    r = Reader(data_package, "results 0.001mm", f"./res/{data_package}")
    r.process_centre("SANPIETRO")


if __name__ == '__main__':
    main()
