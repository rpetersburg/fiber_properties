from FiberProperties import ImageAnalysis, modalNoise, circleArray, plotFFT, showPlot
import numpy as np

height = 2000
width = 2000
mesh_grid = np.meshgrid(np.arange(height), np.arange(width))

fft_list = []
freq_list = []
label_list = []

perfect_image = 65535 * circleArray(mesh_grid, height/2, width/2, 0.8*min(height,width))

for factor in [0.01, 0.1, 1.0]:
    test_image = perfect_image * factor
    baseline_image = np.zeros_like(test_image)
    for i in xrange(10):
        baseline_image += np.random.poisson(test_image) / 10

    baseline_obj = ImageAnalysis(baseline_image,
                                 pixel_size=3.45,
                                 camera='nf')
    fft, freq = modalNoise(baseline_obj,
                           method='fft',
                           output='array',
                           radius_factor=1.0,
                           show_image=False)
    fft_list.append(fft)
    freq_list.append(freq)
    label_list.append('baseline_' + str(factor))

plotFFT(freq_list,
        fft_list,
        labels=label_list,
        title='Modal Noise Baseline Test')
showPlot()