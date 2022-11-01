from MyCacher import *
import multiprocessing
import time

def main():
    start_time = time.time()

    p1 = multiprocessing.Process(target=Cacher, args=("./res/Y_Fields/", "./out/Y_Fields/", 0.1))
    p1.start()

    p1.join()

    end_time = time.time()
    print(f"exec time: {end_time - start_time}")


if __name__ == '__main__':
    main()
