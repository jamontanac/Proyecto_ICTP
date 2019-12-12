import numpy as np
import matplotlib.pylab as plt
import os

names=[i for i in os.listdir(".") if i.endswith(".txt")]
#print(names)
plt.figure(figsize=(5,3*len(names)))
#plt.subplots_adjust(hspace = .001)
aux=1
for i in sorted(names):
	print(i)
	a=np.genfromtxt(i)
	plt.subplot(len(names),1,aux)
	plt.imshow(a,cmap="tab20")
	plt.colorbar()
	plt.clim(0,100)
	aux+=1
plt.savefig("Matrix.png")
plt.close()
# for i in Rango:
# 	a=np.genfromtxt("matrix_"+str(i)+".txt")
# 	plt.subplot(len(Rango),1,i+1)
#
# 	plt.imshow(a)
# 	plt.axis('off')
# 	#plt.colorbar()
# 	plt.clim(0,100)
# plt.savefig("Matrix.png")
# plt.close()
