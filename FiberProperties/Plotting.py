"""Plotting.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The functions in this module are used to plot graphs and images relevant to
the FiberProperties package
"""
import matplotlib.pyplot as plt
plt.rc('font', size=16, family='serif')
plt.rc('figure', figsize=[20, 12.36])
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('lines', lw=4)

#=============================================================================#
#===== General Use Functions =================================================#
#=============================================================================#

def showPlots():
    plt.show()

def savePlot(file):
    plt.savefig(file)

#=============================================================================#
#===== Numpy Array Plotting ==================================================#
#=============================================================================#

def plotHorizontalCrossSection(image_array, row):
    row_int = int(round(row))
    plt.plot(image_array[row_int, :])
    plt.title('Horizontal Cross Section (row = %s)'%row)
    plt.xlabel('Pixel')
  
def plotVerticalCrossSection(image_array, column):
    column_int = int(round(column))
    plt.plot(image_array[:, column_int])
    plt.title('Vertical Cross Section (column = %s)'%column)
    plt.xlabel('Pixel')

def plotCrossSections(image_array, row, column):
    plt.figure()
    plt.subplot(211)
    plotHorizontalCrossSection(image_array, row)
    plt.subplot(212)
    plotVerticalCrossSection(image_array, column)

def plotOverlaidCrossSections(first_array, second_array, row, column):
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

def plotDot(image_array, row, column):
    plotImageArray(image_array)
    plt.scatter(column, row, s=25, color='red')

def plotCrossSectionSums(image_array):
    plt.figure()
    plt.subplot(211)
    plt.plot(sumRows(image_array))
    plt.title('Average for each Column')
    plt.xlabel('Column')
    plt.subplot(212)
    plt.plot(sumColumns(image_array))
    plt.title('Average for each Row')
    plt.xlabel('Row')

def plotImageArray(image_array):
    plt.figure()
    plt.imshow(image_array, cmap='gray')
    plt.colorbar(label='intensity')
    plt.xlabel('x pixel')
    plt.ylabel('y pixel')

def showImageArray(image_array):
    plotImageArray(image_array)
    showPlots()

#=============================================================================#
#===== Fiber Property Plotting ===============================================#
#=============================================================================#

def plotFFT(freq_arrays, fft_arrays, labels=['No label'], title='Power Spectrum', min_wavelength=None, max_wavelength=20.0):
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
