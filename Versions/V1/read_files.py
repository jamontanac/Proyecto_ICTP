import numpy as np
import matplotlib.pylab as plt
plt.figure(figsize=(5,10))
for i in range(0,4):
    a=np.genfromtxt("matrix_"+str(i)+".txt")
    plt.subplot(4,1,i+1)
    plt.imshow(a,cmap="tab20")
    plt.colorbar()
    plt.clim(0,100)
plt.savefig("matrix.png")
plt.close()
