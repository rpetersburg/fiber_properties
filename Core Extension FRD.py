from ImageAnalysis import ImageAnalysis
from Calibration import Calibration
from FiberProperties import FRD
from NumpyArrayHandler import saveArray
import matplotlib.pyplot as plt

def coreExtensionFRD(folder, name):
    """Collects all relevant FRD info for core extension fibers

    Args:
        folder: top level folder where images are contained
        name: name to be used for saving plots

    Returns:
        input_focal_ratio [list]:
        energy_loss [list]:
        output_focal_ratio [list]:
    """
    calibration = Calibration([folder + 'Dark/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                              None,
                              [folder + 'Ambient/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

    images = [folder + 'Output 3.5/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
    diameter = ImageAnalysis(images, calibration, magnification=1, threshold=10000).getFiberDiameter(method='edge', units='microns', show_image=False)
    magnification = diameter / ((4.0 / 3.5) * 25400)
    print 'Magnification:', magnification

    fig1, ax1 = plt.subplots()
    input_focal_ratio = [3.0, 3.5, 4.0, 4.5, 5.0]
    energy_loss = []
    output_focal_ratio = []
    for f in input_focal_ratio:
        print 'F/' + str(f) + ':'
        images = [folder + 'Input ' + str(f) + '/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
        ff_obj = ImageAnalysis(images, calibration, magnification=magnification)
        frd_output = FRD(ff_obj, input_focal_ratio=f, focal_lim=(2.3, 6.0), res=0.1)

        ax1.plot(frd_output[0], frd_output[1], label=str(f))
        energy_loss.append(frd_output[2])
        print 'Energy loss:', frd_output[2], '%'
        output_focal_ratio.append(frd_output[3])
        print 'Output F/#:', frd_output[3]
        saveArray(ff_obj.getImage(), folder + 'Input ' + str(f) + ' Image.tif')

    ax1.set_xlabel('Output f/#')
    ax1.set_ylabel('Encircled Energy')
    ax1.legend(title='Input f/#')
    ax1.grid()
    ax1.set_title('FRD: ' + name)
    fig1.savefig(folder + name + ' FRD.png')

    print

    file = open(folder + 'Info.txt', 'w')
    file.write(str({'focal_ratios': input_focal_ratio,
                    'energy_loss': energy_loss,
                    'output_focal_ratio': output_focal_ratio}))
    file.close()

    return input_focal_ratio, energy_loss, output_focal_ratio

if __name__ == '__main__':
    print 'Reference Fiber:'
    folder = 'FRD Measurements/Core Extension/2016-08-04 Reference Octagonal/'
    ref_info = coreExtensionFRD(folder, 'Reference Fiber')

    print 'Prototype Fiber 1:'
    folder = 'FRD Measurements/Core Extension/2016-08-05 Prototype Core Extension 1/'
    prototype1_info = coreExtensionFRD(folder, 'Prototype Fiber 1')

    print 'Prototype Fiber 2:'
    folder = 'FRD Measurements/Core Extension/2016-08-08 Prototype Core Extension 2/'
    prototype2_info = coreExtensionFRD(folder, 'Prototype Fiber 2')

    input_focal_ratio = ref_info[0]

    plt.figure()
    plt.plot(input_focal_ratio, ref_info[1], label='Reference')    
    plt.plot(input_focal_ratio, prototype1_info[1], label='Prototype 1')
    plt.plot(input_focal_ratio, prototype2_info[1], label='Prototype 2')
    plt.xlabel('Input f/#')
    plt.ylabel('Attenuation [%]')
    plt.grid()
    plt.legend()
    plt.title('Attenuation at Constant F/#')
    plt.savefig('FRD Measurements/Core Extension/Attenuation.png')

    plt.figure()
    plt.plot(input_focal_ratio, ref_info[2], label='Reference')
    plt.plot(input_focal_ratio, prototype1_info[2], label='Prototype 1')
    plt.plot(input_focal_ratio, prototype2_info[2], label='Prototype 2')
    plt.plot(input_focal_ratio, input_focal_ratio, label='Ideal', linestyle='--')
    plt.xlabel('Input f/#')
    plt.ylabel('Output f/#')
    plt.grid()
    plt.legend()
    plt.title('FRD Comparison')
    plt.savefig('FRD Measurements/Core Extension/Input vs Output.png')

