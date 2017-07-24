import numpy as np
from fiber_properties import FiberImage, modal_noise, circle_array, gaussian_array, plot_fft, show_plot

height = 2000
width = 2000
mesh_grid = np.meshgrid(np.arange(height), np.arange(width))
amp = 40000

fft_list = []
freq_list = []
label_list = []

perfect_image = amp * circle_array(mesh_grid, height/2, width/2, 0.8*min(height,width))
# perfect_image = gaussian_array(mesh_grid, height/2, width/2, 0.8*min(height,width), amp, 0).reshape(height, width)

for factor in [0.1, 0.2, 0.4, 0.8, 1.0]:
    test_image = perfect_image * factor
    baseline_image = np.zeros_like(test_image)
    for i in xrange(10):
        baseline_image += np.random.poisson(test_image) / 10

    baseline_obj = FiberImage(baseline_image,
                                 pixel_size=3.45,
                                 camera='nf')
    fft, freq = modal_noise(baseline_obj,
                            method='fft',
                            output='array',
                            radius_factor=1.0,
                            show_image=False)

    fft_list.append(fft)
    freq_list.append(freq)
    label_list.append('baseline_' + str(factor))

plot_fft(freq_list[:-1],
         [fft_list[i] - fft_list[-1] for i in xrange(4)],
         labels=label_list[:-1],
         title='Modal Noise Baseline Test')
show_plot()