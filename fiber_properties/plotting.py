"""Plotting.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The functions in this module are used to plot graphs and images relevant to
the FiberProperties package
"""
import matplotlib.pyplot as plt
from collections import Iterable
import numpy as np
from .numpy_array_handler import sum_rows, sum_columns
from .input_output import create_directory
plt.rc('figure', figsize=[3.39, 3.0])
plt.rc('text', usetex=True)

plt.rc('font', size=10, family='serif', serif=['Computer Modern Roman'])
plt.rc('axes', labelsize=10, linewidth=1)
plt.rc('legend', frameon=True, fontsize=8, labelspacing=0.3, numpoints=1)
plt.rc('lines', linewidth=1)

plt.rc('xtick', labelsize=10)
plt.rc('xtick.major', size=4, width=1)
plt.rc('xtick.minor', visible=True, size=2, width=1)

plt.rc('ytick', labelsize=10)
plt.rc('ytick.major', size=4, width=1)
plt.rc('ytick.minor', visible=True, size=2, width=1)

#=============================================================================#
#===== General Use Functions =================================================#
#=============================================================================#

def show_plots():
    plt.tight_layout()
    plt.show()

def save_plot(file, **kwargs):
    create_directory(file)
    plt.tight_layout()
    plt.savefig(file, **kwargs)

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
    plt.figure(figsize=[6,6])
    plt.subplot(211)
    plot_horizontal_cross_section(image, pixel.y)
    plt.xlabel('')
    plt.subplot(212)
    plot_vertical_cross_section(image, pixel.x)

def plot_overlaid_cross_sections(first_array, second_array, pixel):
    row = int(round(pixel.y))
    column = int(round(pixel.x))
    plt.figure(figsize=[6,6])
    plt.subplot(211)
    plt.plot(first_array[row, :])
    plt.plot(second_array[row, :])
    plt.title('Horizontal Cross Section (row = %s)'%row)
    plt.subplot(212)
    plt.plot(first_array[:, column])
    plt.plot(second_array[:, column])
    plt.title('Vertical Cross Section (column = %s)'%column,)
    plt.xlabel('Pixel')

def plot_dot(image, pixel):
    plot_image(image)
    plt.scatter(pixel.x, pixel.y, s=25, color='red')

def plot_dots(image, *pixels):
    plot_image(image)
    for pixel in pixels:
        plt.scatter(pixel.x, pixe.y, s-25, color='red')

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
    plt.ylabel('x drift ($\mu m$)', fontsize='small')
    plt.title(r'$\sigma_{%s}= %.3f um, SG_{max} = %d$' % (cam, sigma, max_sg))
    plt.xlim(min(data.time), max(data.time))

    plt.subplot(312)
    plt.plot(data.time, data.y_diff)
    plt.ylabel('y drift ($\mu m$)', fontsize='small')
    plt.xlim(min(data.time), max(data.time))

    plt.subplot(313)
    plt.plot(data.time, data.diameter)
    plt.xlabel('time [min]')
    plt.ylabel('diameter ($\mu m$)', fontsize='small')
    plt.xlim(min(data.time), max(data.time))

def plot_stability_binned(data, cam, bin_size):
    plot_stability(data, cam)

    half_bin = bin_size / 2
    time = data.time[half_bin:len(data.time)-half_bin+1:bin_size]
    x_diff = []
    y_diff = []
    x_std = []
    y_std = []
    diameter = []
    for i in xrange(half_bin, len(data.time)-half_bin+1, bin_size):
        x_data = np.array(data.x_diff[i-half_bin:i+half_bin])
        y_data = np.array(data.y_diff[i-half_bin:i+half_bin])
        x_diff.append(x_data.mean())
        y_diff.append(y_data.mean())
        x_std.append(x_data.std())
        y_std.append(y_data.std())
        diameter.append(np.array(data.diameter[i-half_bin:i+half_bin]).mean())

    # sigma = np.sqrt(np.std(x_diff)**2 + np.std(y_diff)**2)
    sigma = np.sqrt(np.mean(x_std)**2 + np.mean(y_std)**2)
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

def plot_modal_noise(modal_noise, plot_type='bar', labels=[''], bar_labels=[], method='filter', errors=[[]]):
    plt.figure()
    colors = ['C3', 'C1', 'C8', 'C2', 'C0', 'C4', 'C6', 'C5', 'C9', 'C7']
    num_tests = len(labels)

    if plot_type is 'bar':
        if not bar_labels:
            raise RuntimeError('Please include labels for modal noise bar plot')

        bar_width = 0.8 / num_tests
        indexes = np.arange(len(bar_labels))

        plt.grid(which='major', axis='y', zorder=0)

        for i, (mn, label, color) in enumerate(zip(modal_noise, labels, colors)):
            plt.bar(indexes+0.1+i*bar_width,
                    mn, bar_width, label=label,
                    color=color, edgecolor='none',
                    zorder=3)
            if errors and errors[0]:
                plt.errorbar(indexes+0.1+i*bar_width, mn, yerr=errors[i],
                             ecolor='k', zorder=5, fmt='none')


        plt.xticks(indexes+0.5, bar_labels, rotation=30, ha='right')
        plt.tick_params(axis='x', which='both', bottom='off', top='off')

    elif plot_type is 'line':
        for mn, label, color in zip(modal_noise, labels, colors):
            plt.plot(range(1, len(mn)+1), mn, label=label, color=color, marker='o')
            if errors:
                plt.errorbar(range(1, len(mn)+1), mn, yerr=errors[i], ecolor=color, fmt='none')
        plt.xlabel('Integration Time [s]')

    else:
        raise RuntimeError('Please input valid plot_type')

    if num_tests > 1:
        plt.legend(loc='best', frameon=True)

    y_label = 'Signal to Noise Ratio'
    if method == 'contrast':
        y_label = 'Michelson Contrast'
    elif method == 'gradient':
        y_label = 'Gradient NSR'
    elif method == 'gini':
        y_label = 'Gini Coefficient'
    elif method == 'entropy':
        y_label = 'Hartley Entropy'
    plt.ylabel(y_label)
    plt.margins(0.1)
    plt.ylim(ymin=0.0)

def plot_fft(fft_info, labels=[],
             min_wavelength=None, max_wavelength=100.0):
    plt.figure()
    power = []

    if isinstance(fft_info, Iterable):
        if not labels:
            raise RuntimeError('Please label power spectra')
        for info, label in zip(fft_info, labels):
            _plot_fft(info, label)
            power.append(info.power)
    else:
        info = fft_info
        _plot_fft(info, labels)
        power.append(info.power)

    if min_wavelength is None:
        min_wavelength = 2.0/info.freq.max()
    min_power = np.array(power).min()

    plt.xlim(min_wavelength, max_wavelength)
    plt.ylim(ymin=min_power*0.9)
    plt.xscale('log', subsx=[2,3,4,5,6,7,8,9])
    plt.yscale('log', subsy=[2,3,4,5,6,7,8,9])
    plt.ylabel('normalized power')
    plt.xlabel('speckle diameter ($\mu m$)')
    # plt.grid(which='both', axis='x')
    if labels:
        plt.legend(loc=2, frameon=False)

def _plot_fft(fft_info, label):
    fft_info.freq[0] += fft_info.freq[1] / 10.0
    wavelength_array = 1.0 / fft_info.freq
    if label:
        plt.plot(wavelength_array, fft_info.power, label=label)
    else:
        plt.plot(wavelength_array, fft_info.power)

def plot_scrambling_gain_input_output(scrambling_output):
    input_x = scrambling_output.in_x
    input_y = scrambling_output.in_y
    output_x = scrambling_output.out_x
    output_y = scrambling_output.out_y
    plt.figure()
    plt.subplot(221)
    plt.scatter(input_x, output_x)
    plt.xlabel('$x/D_{in}$')
    plt.ylabel('$x/D_{out}$')
    plt.subplot(222)
    plt.scatter(input_y, output_x)
    plt.xlabel('$y/D_{in}$')
    plt.ylabel('$x/D_{out}$')
    plt.subplot(223)
    plt.scatter(input_x, output_y)
    plt.xlabel('$x/D_{in}$')
    plt.ylabel('$y/D_{out}$')
    plt.subplot(224)
    plt.scatter(input_y, output_y)
    plt.xlabel('$y/D_{in}$')
    plt.ylabel('$y/D_{out}$')

def plot_scrambling_gain(scrambling_output):
    input_dist = scrambling_output.in_d
    output_dist = scrambling_output.out_d
    scrambling_gain = scrambling_output.scrambling_gain
    plt.figure()
    plt.subplot(211)
    plt.scatter(input_dist, output_dist)
    plt.xlabel('$\Delta d_{in}/D_{in}')
    plt.ylabel('$\Delta d_{out}/D_{out}')

    plt.subplot(212)
    plt.scatter(input_dist, scrambling_gain)
    plt.xlabel('$\Delta d_{in}/D_{in}')
    plt.ylabel('scrambling gain')

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
    plt.xlabel('output f/#')
    plt.ylabel('encircled energy')
    plt.ylim(ymax=1)
    plt.grid()
    plt.legend(loc=3, title='input f/#')

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
    plt.xlabel('input f/#')
    plt.ylabel('energy loss (\
        %)')
    plt.grid()
    plt.legend(loc=2)

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
    plt.xlabel('input f/#')
    plt.ylabel('output f/#')
    plt.grid()
    plt.legend(loc=2)

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
                plt.xlabel('output f/#')
                plt.ylabel('encircled energy')
                plt.ylim(ymax=1)
                plt.grid()
                plt.legend(title='input f/' + str(f), loc=3)
