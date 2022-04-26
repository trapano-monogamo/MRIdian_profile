from MyCacher import *

def main():
   # c = Cacher("./res/", "./output tests/original_not_filtered/", TestSettings(0, moving_average, [3]))
   # c = Cacher("./res/", "./output tests/dose_filtered - spline/", TestSettings(1, "spline"))
   # c = Cacher("./res/", "./output tests/first_derivative_filtered/", TestSettings(2, moving_average))
   # c = Cacher("./res/", "./output tests/dose_first_derivative_filtered/", TestSettings(3, "moving_average"))
   # c = Cacher("./res/", "./output tests/second_derivative_filtered/", TestSettings(4, moving_average, [3]))
   # c = Cacher("./res/", "./output tests/both_derivatives_filtered/", TestSettings(5, moving_average, [3]))
   # c = Cacher("./res/", "./output tests/all_filtered/", TestSettings(6, moving_average, [3]))

   c = Cacher("./res/", "./out/sanpietro 0.1mm - MA 13 - dose_filtered/", TestSettings(1, "moving_average"))

if __name__ == '__main__':
   main()
