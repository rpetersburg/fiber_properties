from ImageAnalysis import ImageAnalysis
from Calibration import Calibration
from FiberProperties import FRD
import matplotlib.pyplot as plt

def coreExtensionFRD(folder, name):
    """Collects all relevant FRD info for core extension fibers

    Args:
        folder: top level folder where images are contained
        name: name to be used for saving plots

    Returns:
        input_f_number [list]:
        energy_loss [list]:
        output_f_number [list]:
    """
    calibration = Calibration([folder + 'Dark/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                              None,
                              [folder + 'Ambient/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

    images = [folder + 'Output 3.5/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
    diameter = ImageAnalysis(images, calibration, magnification=1, threshold=10000).getFiberDiameter(method='edge', units='pixels', show_image=False)
    magnification = (diameter * 3.45) / ((4.0 / 3.5) * 25400)
    print 'Magnification:', magnification

    plt.figure(1)
    input_f_number = [3.0, 3.5, 4.0, 4.5, 5.0]
    energy_loss = []
    output_f_number = []
    for f in input_f_number:
        print 'F/' + str(f) + ':'
        images = [folder + 'Input ' + str(f) + '/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
        ff_obj = ImageAnalysis(images, calibration, magnification=magnification)
        frd_output = FRD(ff_obj, input_f_number=f, f_lim=(2.4, 10.0), res=0.1)

        plt.plot(frd_output[0], frd_output[1], label=str(f))
        energy_loss.append(frd_output[2])
        print 'Energy loss:', frd_output[2], '%'
        output_f_number.append(frd_output[3])
        print 'Output F/#:', frd_output[3]

    plt.xlabel('Output f/#')
    plt.ylabel('Encircled Energy')
    plt.legend(title='Input f/#')
    plt.grid()
    plt.title('FRD: ' + name)
    plt.savefig(folder + name + ' FRD.png')
    plt.clf()
    plt.cla()
    plt.close()

    print

    return input_f_number, energy_loss, output_f_number

if __name__ == '__main__':
    print 'Prototype Fiber 1:'
    folder = 'FRD Measurements/Core Extension/2016-08-05 Prototype Core Extension 1/'
    prototype1_info = coreExtensionFRD(folder, 'Prototype Fiber 1')

    print 'Reference Fiber:'
    folder = 'FRD Measurements/Core Extension/2016-08-04 Reference Octagonal/'
    ref_info = coreExtensionFRD(folder, 'Reference Fiber')

    input_f_number = ref_info[0]

    plt.figure(2)
    plt.plot(input_f_number, prototype1_info[1], label='Prototype 1')
    plt.plot(input_f_number, ref_info[1], label='Reference')
    plt.xlabel('Input f/#')
    plt.ylabel('Attenuation [%]')
    plt.grid()
    plt.legend()
    plt.title('Attenuation at Constant F/#')
    plt.savefig('FRD Measurements/Core Extension/Attenuation.png')

    plt.figure(3)
    plt.plot(input_f_number, prototype1_info[2], label='Prototype 1')
    plt.plot(input_f_number, ref_info[2], label='Reference')
    plt.plot(input_f_number, input_f_number, label='Ideal', linestyle='-')
    plt.xlabel('Input f/#')
    plt.ylabel('Output f/#')
    plt.grid()
    plt.legend()
    plt.title('FRD Comparison')
    plt.savefig('FRD Measurements/Core Extension/Input vs Output.png')

