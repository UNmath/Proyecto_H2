import numpy as np
import matplotlib.pyplot as plt

datos = np.loadtxt('datos.txt', delimiter = ',')

#y = (1/datos[:,1]) * np.sin(datos[:,0])

#x = (1/datos[:,1]) * np.cos(datos[:,0])

x = datos[:,2]
y = datos[:,3]

plt.figure()
plt.plot(x,y)
plt.axis('equal')
plt.title('Trayectoria fotón en agujero de kerr en el ecuador')
plt.savefig('figura.png')