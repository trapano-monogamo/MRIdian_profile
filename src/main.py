from MyCacher import *

def main():
   c = Cacher("./res/", "./out/")

   # x = np.linspace(0, len(test_scan.derivative), len(test_scan.derivative))
   # y = test_scan.derivative

   # peak_max = max(test_scan.derivative)
   # pmaxi = test_scan.derivative.index(peak_max)
   # peak_min = min(test_scan.derivative)
   # pmini = test_scan.derivative.index(peak_min)

   # ext_data = []
   # for i in range(len(test_scan.data) - 1):
   #     ext_data.append(test_scan.data[i][1])

   # plt.plot(x, ext_data, c = "blue")
   # plt.plot(x, y, c = "red")
   # plt.scatter([pmaxi, pmini], [peak_max, peak_min], c= "green")
   # plt.show()

   # MAKE ALL PROFILE AUTOMATICALLY PRODUCE PLOTS
   # test_scan.find_inflection_points()
   # test_scan.output_plot()

if __name__ == '__main__':
   main()
