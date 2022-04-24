from MyCacher import *

def main():
   # c = Cacher("./res/", "./output tests/original_not_filtered/", TestSettings(moving_average, 3, 0))
   # c = Cacher("./res/", "./output tests/dose_filtered/", TestSettings(moving_average, 3, 1))
   # c = Cacher("./res/", "./output tests/first_derivative_filtered/", TestSettings(moving_average, 3, 2))
   # c = Cacher("./res/", "./output tests/second_derivative_filtered/", TestSettings(moving_average, 3, 3))
   # c = Cacher("./res/", "./output tests/both_derivatives_filtered/", TestSettings(moving_average, 3, 4))
   c = Cacher("./res/", "./output tests/all_filtered/", TestSettings(moving_average, 3, 5))

   # c = Cacher("./res/", "./out/", TestSettings(moving_average, 3, 4))

if __name__ == '__main__':
   main()
