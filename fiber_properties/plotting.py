"""Plotting.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The functions in this module are used to plot graphs and images relevant to
the FiberProperties package
"""
import matplotlib.pyplot as plt
from fiber_properties.numpy_array_handler import sum_rows, sum_columns
plt.rc('font', size=16, family='serif')
plt.rc('figure', figsize=[20, 12.36])
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('lines', lw=4)

#=============================================================================#
#===== General Use Functions =================================================#
#=============================================================================#

def show_plots():
    plt.show()

def save_plot(file):
    plt.savefig(file)

#=============================================================================#
#===== Numpy Array Plotting ==================================================#
#=============================================================================#

def plot_horizontal_cross_section(image_array, row):
    row_int = int(round(row))
    plt.plot(image_array[row_int, :])
    plt.title('Horizontal Cross Section (row = %s)'%row)
    plt.xlabel('Pixel')

def plot_vertical_cross_section(image_array, column):
    column_int = int(round(column))
    plt.plot(image_array[:, column_int])
    plt.title('Vertical Cross Section (column = %s)'%column)
    plt.xlabel('Pixel')

def plot_cross_sections(image_array, row, column):
    plt.figure()
    plt.subplot(211)
    plot_horizontal_cross_section(image_array, row)
    plt.subplot(212)
    plot_vertical_cross_section(image_array, column)

def plot_overlaid_cross_sections(first_array, second_array, row, column):
    row = int(round(row))
    column = int(round(column))
    plt.figure()
    plt.subplot(211)
    plt.plot(first_array[row, :])
    plt.plot(second_array[row, :])
    plt.title('Horizontal Cross Section (row = %s)'%row)
    plt.xlabel('Pixel')
    plt.subplot(212)
    plt.plot(first_array[:, column])
    plt.plot(second_array[:, column])
    plt.title('Vertical Cross Section (column = %s)'%column,)
    plt.xlabel('Pixel')

def plot_dot(image_array, row, column):
    plot_image_array(image_array)
    plt.scatter(column, row, s=25, color='red')

def plot_cross_section_sums(image_array):
    plt.figure()
    plt.subplot(211)
    plt.plot(sum_rows(image_array))
    plt.title('Average for each Column')
    plt.xlabel('Column')
    plt.subplot(212)
    plt.plot(sum_columns(image_array))
    plt.title('Average for each Row')
    plt.xlabel('Row')

def plot_image_array(image_array):
    plt.figure()
    plt.imshow(image_array, cmap='gray')
    plt.colorbar(label='intensity')
    plt.xlabel('x pixel')
    plt.ylabel('y pixel')

def show_image_array(image_array):
    plot_image_array(image_array)
    show_plots()

#=============================================================================#
#===== Fiber Property Plotting ===============================================#
#=============================================================================#

def plot_fft(freq_arrays, fft_arrays, labels=['No label'],
             title='Power Spectrum', min_wavelength=None, max_wavelength=20.0):
    plt.figure()
    wavelength_arrays = []
    for i in xrange(len(freq_arrays)):
        wavelength_arrays.append(1.0/freq_arrays[i])
        plt.plot(wavelength_arrays[i], fft_arrays[i], label=labels[i])
    if min_wavelength is None:
        min_wavelength = 2.0/freq_arrays[0].max()
    plt.xlim(min_wavelength, max_wavelength)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('Normalized Power')
    plt.xlabel('Speckle Diameter [um]')
    plt.title(title)
    plt.legend()

def plot_scrambling_gain_input_output(scrambling_output, camera=''):
    input_x = scrambling_output[0]
    input_y = scrambling_output[1]
    output_x = scrambling_output[2]
    output_y = scrambling_output[3]
    plt.figure()
    plt.subplot(221)
    plt.scatter(input_x, output_x)
    plt.xlabel('Input X [Fiber Diameter]')
    plt.ylabel('Output X [Fiber Diameter]')
    plt.subplot(222)
    plt.scatter(input_y, output_x)
    plt.xlabel('Input Y [Fiber Diameter]')
    plt.ylabel('Output X [Fiber Diameter]')
    plt.subplot(223)
    plt.scatter(input_x, output_y)
    plt.xlabel('Input X [Fiber Diameter]')
    plt.ylabel('Output Y [Fiber Diameter]')
    plt.subplot(224)
    plt.scatter(input_y, output_y)
    plt.xlabel('Input Y [Fiber Diameter]')
    plt.ylabel('Output Y [Fiber Diameter]')
    plt.suptitle(title + ' Centroid Shift')

def plot_scrambling_gain(scrambling_output, camera=''):
    input_dist = scrambling_output[5]
    output_dist = scrambling_output[6]
    plt.figure()
    plt.scatter(input_dist, output_dist)
    plt.xlabel('Input Delta [Fiber Diameter]')
    plt.ylabel('Output Delta [Fiber Diameter]')
    plt.title('NF Scrambling Gains')
