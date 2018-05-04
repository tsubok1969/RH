import math
import Hau
import numpy as np
import matplotlib.pyplot as plt

beta, theta = map(float, input("Input the upstream beta and theta\n").split())

ax2 = np.arange(100)*0.1
ax1 = Hau.ratio(beta, theta, ax2)

plt.xlim(0.,10.)
plt.ylim(0.,10.)
plt.plot(ax2, ax1)
plt.plot(ax2, ax2)
plt.show()
