import time
from scipy.signal import medfilt2d, medfilt, order_filter
import numpy as np
from fiber_properties import FiberImage, filter_image, plot_image, show_plots

# a = 2**16 * np.random.rand(2200,2200)
a = FiberImage('../data/modal_noise/Kris_data/circular_100um/unagitated/ff_corrected.fit').get_image()
plot_image(a)

start = time.time()
b = filter_image(a, 101, quick=True)
end = time.time()
print(end-start)
plot_image(b[100:200,100:200])

start = time.time()
c = filter_image(a, 101, quick=False)
end = time.time()
print(end-start)
plot_image(c[100:200,100:200])

show_plots()
