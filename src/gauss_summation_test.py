from MyUtils import *
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import random
from scipy.optimize import curve_fit, minimize
from scipy.stats import rv_continuous
from lmfit import Model

matplotlib.use('Svg')

a = 0.9
b = -0.8
c = 6.0

initial_parameters = [0.1, -0.1, 0.2]

x = np.linspace(-2,2,200)
y = [skew_normal(z,a,b,c) for z in x]
y1 = [n + random.uniform(-0.1, 0.1) for n in y]

fargs, pcov = curve_fit(skew_normal, x, y1, initial_parameters)
y2 = [skew_normal(z, *fargs) for z in x]

# model = Model(skew_normal)
# params = model.make_params(a=0.1, b=-0.1, c=0.1)
# y2 = model.eval(params, x=x)

# params = minimize(skew_normal, np.array([0.1,-0.1,0.1]), args=(x,y1,x), method='Nelder-Mead')
# y2 = [skew_normal(z, *(params.x)) for z in x]

# plt.plot(x,y)
plt.plot(x,y1)
plt.plot(x,y2)
plt.savefig("./src/test.png")
