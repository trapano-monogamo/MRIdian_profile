from MyCacher import *

def main():
   Cacher(
      "./res/MR_LINAC/",
      f"./out/MR_LINAC/chi2 0.001mm/",
      0.001)
   Cacher(
      "./res/TRUEBEAM/",
      "./out/TRUEBEAM/chi2 0.001mm/",
      0.001)

if __name__ == '__main__':
   main()
