import numpy as np
import matplotlib.pyplot as plt
    
var_array = np.arange(1000) / 10.0
x, y = np.meshgrid(var_array, var_array)

array = np.sin(x) + np.sin(y)
fft = np.abs(np.fft.fft2(array, norm='ortho'))**2

plt.figure(1)
plt.imshow(array, cmap='gray')
plt.colorbar(label='intensity')

plt.figure(2)
plt.imshow(fft, cmap='gray')
plt.colorbar(label='intensity')

plt.show()