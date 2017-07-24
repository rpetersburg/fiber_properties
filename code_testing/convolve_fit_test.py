from fiber_properties import FiberImage, circle_array, plot_image, show_plots
from scipy.ndimage.filters import convolve
import numpy as np

im = FiberImage('../data/modal_noise/amp_freq_200um/baseline/nf_corrected.fit')
image = im.get_image()
diameter = im.get_fiber_diameter(method='edge')
edge = diameter + 10

mask = circle_array(np.meshgrid(np.arange(edge), np.arange(edge)), edge/2.0, edge/2.0, diameter/2.0, res=16)

a = convolve(image, mask)

plot_image(a)
show_plots()