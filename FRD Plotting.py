import matplotlib.pyplot as plt
import numpy as np

def encircledEnergy(input_list):
    return np.array(input_list) / max(input_list)

out = [2, 2.5, 3, 3.5, 4, 4.5, 5]
in30 = [43.5, 43.3, 43.3, 43.1, 42.8, 41.7, 40.1]
in35 = [43.8, 43.6, 43.5, 43.1, 42.8, 42.0, 40.1]
in40 = [44.3, 44.0, 43.8, 43.5, 43.2, 42.3, 40.6]
in45 = [43.2, 43.1, 43.0, 42.7, 42.4, 41.6, 39.8]
in50 = [41.0, 40.9, 40.8, 40.6, 40.4, 40.0, 38.4]

in30 = encircledEnergy(in30)
in35 = encircledEnergy(in35)
in40 = encircledEnergy(in40)
in45 = encircledEnergy(in45)
in50 = encircledEnergy(in50)

plt.figure(1)
plt.plot(out, in30, label='3.0')
plt.plot(out, in35, label='3.5')
plt.plot(out, in40, label='4.0')
plt.plot(out, in45, label='4.5')
plt.plot(out, in50, label='5.0')
plt.legend(title='Input F/#')
plt.xlabel('Output F/#')
plt.ylabel('Encircled Energy')
plt.show()