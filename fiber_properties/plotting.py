"""Plotting.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The functions in this module are used to plot graphs and images relevant to
the FiberProperties package
"""
import matplotlib.pyplot as plt
from collections import Iterable
from .numpy_array_handler import sum_rows, sum_columns
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

def plot_horizontal_cross_section(image, row):
    row_int = int(round(row))
    plt.plot(image[row_int, :])
    plt.title('Horizontal Cross Section (row = %s)'%row)
    plt.xlabel('Pixel')

def plot_vertical_cross_section(image, column):
    column_int = int(round(column))
    plt.plot(image[:, column_int])
    plt.title('Vertical Cross Section (column = %s)'%column)
    plt.xlabel('Pixel')

def plot_cross_sections(image, row, column):
    plt.figure()
    plt.subplot(211)
    plot_horizontal_cross_section(image, row)
    plt.subplot(212)
    plot_vertical_cross_section(image, column)

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

def plot_dot(image, row, column):
    plot_image(image)
    plt.scatter(column, row, s=25, color='red')

def plot_cross_section_sums(image):
    plt.figure()
    plt.subplot(211)
    plt.plot(sum_rows(image))
    plt.title('Average for each Column')
    plt.xlabel('Column')
    plt.subplot(212)
    plt.plot(sum_columns(image))
    plt.title('Average for each Row')
    plt.xlabel('Row')

def plot_image(image):
    plt.figure()
    plt.imshow(image, cmap='gray')
    plt.colorbar(label='intensity')
    plt.xlabel('x pixel')
    plt.ylabel('y pixel')

def show_image(image):
    plot_image(image)
    show_plots()

#=============================================================================#
#===== Fiber Property Plotting ===============================================#
#=============================================================================#

def plot_fft(fft_info, labels=[], title='Power Spectrum',
             min_wavelength=None, max_wavelength=20.0):
    plt.figure()
    wavelength_arrays = []
    if isinstance(fft_info, Iterable):
        for i, info in enumerate(fft_info):
            wavelength_array = 1.0 / info.freq
            if labels:
                plt.plot(wavelength_array, info.power, label=labels[i])
            else:
                raise RuntimeError('Please label power spectra')
        if min_wavelength is None:
            min_wavelength = 2.0/fft_info[0].freq.max()
    else:
        wavelength_array = 1.0 / fft_info.freq
        if labels:
            plt.plot(wavelength_array, fft_info.power, label=labels)
        else:
            plt.plot(wavelength_array, fft_info.power)
        if min_wavelength is None:
            min_wavelength = 2.0/fft_info.freq.max()
    plt.xlim(min_wavelength, max_wavelength)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('Normalized Power')
    plt.xlabel('Speckle Diameter [um]')
    plt.title(title)
    if labels:
        plt.legend(loc='best')

def plot_scrambling_gain_input_output(scrambling_output, title):
    input_x = scrambling_output.in_x
    input_y = scrambling_output.in_y
    output_x = scrambling_output.out_x
    output_y = scrambling_output.out_y
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
    plt.suptitle(title)

def plot_scrambling_gain(scrambling_output, title):
    input_dist = scrambling_output.in_d
    output_dist = scrambling_output.out_d
    scrambling_gain = scrambling_output.scrambling_gain
    plt.figure()
    plt.subplot(211)
    plt.scatter(input_dist, output_dist)
    plt.xlabel('Input Delta [Fiber Diameter]')
    plt.ylabel('Output Delta [Fiber Diameter]')

    plt.subplot(212)
    plt.scatter(input_dist, scrambling_gain)
    plt.xlabel('Input Delta [Fiber Diameter]')
    plt.ylabel('Scrambling Gain')

    plt.suptitle(title)
