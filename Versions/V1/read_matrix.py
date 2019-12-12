import numpy as np
import matplotlib.pylab as plt
a=np.genfromtxt("matrix.txt")
plt.imshow(a,cmap="tab20")
plt.colorbar()
plt.clim(0,100)
plt.show()
