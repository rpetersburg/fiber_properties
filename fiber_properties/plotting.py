"""Plotting.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The functions in this module are used to plot graphs and images relevant to
the FiberProperties package
"""
import matplotlib.pyplot as plt
from collections import Iterable
import numpy as np
from .numpy_array_handler import sum_rows, sum_columns
plt.rc('figure', figsize=[9,6], titlesize='large', titleweight='bold')

plt.rc('font', size=16, family='serif')
plt.rc('axes', labelsize='large', labelweight='bold')
plt.rc('legend', frameon=True, fontsize='large')
plt.rc('lines', lw=4)

plt.rc('xtick', labelsize='large')
plt.rc('xtick.major', size=16)
plt.rc('xtick.minor', visible=True, size=8)

plt.rc('ytick', labelsize='large')
plt.rc('ytick.major', size=16)
plt.rc('ytick.minor', visible=True, size=8)

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

def plot_cross_sections(image, pixel):
    plt.figure(figsize=[20, 12.36])
    plt.subplot(211)
    plot_horizontal_cross_section(image, pixel.y)
    plt.subplot(212)
    plot_vertical_cross_section(image, pixel.x)

def plot_overlaid_cross_sections(first_array, second_array, pixel):
    row = int(round(pixel.y))
    column = int(round(pixel.x))
    plt.figure(figsize=[20, 12.36])
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

def plot_dot(image, pixel):
    plot_image(image)
    plt.scatter(pixel.x, pixel.y, s=25, color='red')

def plot_cross_section_sums(image):
    plt.figure(figsize=[20, 12.36])
    plt.subplot(211)
    plt.plot(sum_rows(image))
    plt.title('Average for each Column')
    plt.xlabel('Column')
    plt.subplot(212)
    plt.plot(sum_columns(image))
    plt.title('Average for each Row')
    plt.xlabel('Row')

def plot_image(image):
    plt.figure(figsize=[20, 12.36])
    plt.imshow(image, cmap='inferno')
    plt.colorbar(label='intensity')
    plt.xlabel('x pixel')
    plt.ylabel('y pixel')

def show_image(image):
    plot_image(image)
    show_plots()

#=============================================================================#
#===== FCS Stability Plotting ================================================#
#=============================================================================#

def plot_stability(data, cam):
    sigma = np.sqrt(np.std(data.x_diff)**2 + np.std(data.y_diff)**2)
    max_sg = np.median(data.diameter) / sigma

    plt.figure(figsize=[15,10])

    plt.subplot(311)
    plt.plot(data.time, data.x_diff)
    plt.ylabel('x drift [um]', fontsize='small')
    plt.title(r'$\sigma_{%s}= %.3f um, SG_{max} = %d$' % (cam, sigma, max_sg))
    plt.xlim(min(data.time), max(data.time))

    plt.subplot(312)
    plt.plot(data.time, data.y_diff)
    plt.ylabel('y drift [um]', fontsize='small')
    plt.xlim(min(data.time), max(data.time))

    plt.subplot(313)
    plt.plot(data.time, data.diameter)
    plt.xlabel('time [min]')
    plt.ylabel('diameter [um]', fontsize='small')
    plt.xlim(min(data.time), max(data.time))

def plot_stability_binned(data, cam, bin_size):
    plot_stability(data, cam)

    half_bin = bin_size / 2
    time = data.time[half_bin:len(data.time)-half_bin]
    x_diff = []
    y_diff = []
    diameter = []
    for i in xrange(half_bin, len(data.time)-half_bin):
        x_diff.append(np.array(data.x_diff[i-half_bin:i+half_bin]).mean())
        y_diff.append(np.array(data.y_diff[i-half_bin:i+half_bin]).mean())
        diameter.append(np.array(data.diameter[i-half_bin:i+half_bin]).mean())

    sigma = np.sqrt(np.std(x_diff)**2 + np.std(y_diff)**2)
    max_sg = np.median(diameter) / sigma
    print cam, 'max SG:', max_sg

    plt.subplot(311)
    plt.plot(time, x_diff, c='red')
    plt.title(r'$\sigma_{%s}= %.3f um, SG_{max} = %d$' % (cam, sigma, max_sg))

    plt.subplot(312)
    plt.plot(time, y_diff, c='red')

    plt.subplot(313)
    plt.plot(time, diameter, c='red')

#=============================================================================#
#===== Fiber Property Plotting ===============================================#
#=============================================================================#

def plot_fft(fft_info, labels=[],
             min_wavelength=None, max_wavelength=100.0):
    plt.figure()
    wavelength_arrays = []
    if isinstance(fft_info, Iterable):
        for i, info in enumerate(fft_info):
            info.freq[0] += info.freq[1] / 2.0
            wavelength_array = 1.0 / info.freq
            if labels:
                plt.plot(wavelength_array, info.power, label=labels[i])
            else:
                raise RuntimeError('Please label power spectra')
        if min_wavelength is None:
            min_wavelength = 2.0/fft_info[0].freq.max()
    else:
        fft_info.freq[0] += fft_info.freq[1] / 2.0
        wavelength_array = 1.0 / fft_info.freq
        if labels:
            plt.plot(wavelength_array, fft_info.power, label=labels)
        else:
            plt.plot(wavelength_array, fft_info.power)
        if min_wavelength is None:
            min_wavelength = 2.0/fft_info.freq.max()
    plt.xlim(min_wavelength, max_wavelength)
    plt.xscale('log', subsx=[2,3,4,5,6,7,8,9])
    plt.yscale('log', subsy=[2,3,4,5,6,7,8,9])
    plt.ylabel('Normalized Power')
    plt.xlabel('Speckle Diameter [um]')
    # plt.grid(which='both', axis='x')
    if labels:
        plt.legend(loc=2, frameon=False)
    plt.tight_layout()

def plot_scrambling_gain_input_output(scrambling_output):
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
    plt.tight_layout()

def plot_scrambling_gain(scrambling_output):
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

    plt.tight_layout()

def plot_frd_encircled_energy(frd_output):
    frd_info = frd_output[0]
    magn = frd_output[1]
    magn_list = frd_output[2]
    magn_error = frd_output[3]

    plt.figure()
    for i, f in enumerate(frd_info.input_fnum):
        plt.plot(frd_info.encircled_energy_fnum[i],
                 frd_info.encircled_energy[i],
                 label=str(f),
                 linewidth=2)
    plt.xlabel('Output f/#')
    plt.ylabel('Encircled Energy')
    plt.ylim(ymax=1)
    plt.grid()
    plt.legend(loc=3, title='Input f/#')
    plt.tight_layout()

def plot_frd_energy_loss(frd_outputs, labels):
    plt.figure()

    for i, output in enumerate(frd_outputs):
        frd_info = output[0]
        magn = output[1]
        magn_list = output[2]
        magn_error = output[3]

        plt.errorbar(frd_info.input_fnum,
                     frd_info.energy_loss,
                     xerr=magn_error*np.array(frd_info.input_fnum),
                     label=labels[i],
                     linewidth=2)
    plt.xlabel('Input f/#')
    plt.ylabel('Energy Loss [%]')
    plt.grid()
    plt.legend(loc=2)
    plt.tight_layout()

def plot_frd_input_output(frd_outputs, labels, ideal=True):
    plt.figure()

    for i, output in enumerate(frd_outputs):
        frd_info = output[0]
        magn = output[1]
        magn_list = output[2]
        magn_error = output[3]

        plt.errorbar(frd_info.input_fnum,
                     frd_info.output_fnum,
                     xerr=magn_error*np.array(frd_info.input_fnum),
                     yerr=magn_error*np.array(frd_info.input_fnum),
                     label=labels[i],
                     linewidth=2)

    if ideal:
        plt.plot(frd_info.input_fnum, frd_info.input_fnum,
                 label='Ideal', linestyle='--', color='black')
    plt.xlabel('Input f/#')
    plt.ylabel('Output f/#')
    plt.grid()
    plt.legend(loc=2)
    plt.tight_layout()

def plot_frd_encircled_energy_comparison(frd_outputs, labels):
    plt.figure(figsize=[18,18])

    for i, output in enumerate(frd_outputs):
        frd_info = output[0]
        magn = output[1]
        magn_list = output[2]
        magn_error = output[3]

        for j, f in enumerate([2.5, 3.0, 3.5, 4.0, 4.5, 5.0]):
            if f in frd_info.input_fnum:
                plt.subplot(3, 2, j+1)
                index = frd_info.input_fnum.index(f)
                plt.plot(frd_info.encircled_energy_fnum[index],
                         frd_info.encircled_energy[index],
                         label=labels[i],
                         linewidth=2)
                plt.xlabel('Output f/#')
                plt.ylabel('Encircled Energy')
                plt.ylim(ymax=1)
                plt.grid()
                plt.legend(title='Input f/' + str(f), loc=3)
    plt.tight_layout()
