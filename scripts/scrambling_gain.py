from FiberProperties import loadImageObject, scramblingGain
import matplotlib.pyplot as plt

if __name__ == '__main__':

    NEW_DATA = False
    folder = '../data/scrambling/2016-08-05 Prototype Core Extension 1/'

    if NEW_DATA:
        from FiberProperties import ImageAnalysis, Calibration

        in_calibration = Calibration(dark=[folder + 'Dark/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                                     ambient=[folder + 'Ambient/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

        nf_calibration = Calibration(dark=[folder + 'Dark/nf_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                                     ambient=[folder + 'Ambient/nf_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

        ff_calibration = Calibration(dark=[folder + 'Dark/ff_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                                     ambient=[folder + 'Ambient/ff_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

        for shift in ['00', '05', '10', '15', '20', '25', '30']:
            print 'Initializing Shift ' + shift
            in_images = [folder + 'Shift_' + shift + '/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
            nf_images = [folder + 'Shift_' + shift + '/nf_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
            ff_images = [folder + 'Shift_' + shift + '/ff_' + str(i).zfill(3) + '.fit' for i in xrange(10)]

            ImageAnalysis(in_images, calibration=in_calibration, camera='in').save()
            ImageAnalysis(nf_images, calibration=nf_calibration, camera='nf').save()
            ImageAnalysis(ff_images, calibration=ff_calibration, camera='ff').save()

    shifts = ['00', '05', '10', '15', '20', '25', '30']
    in_objs = [folder + 'Shift_' + shift + '/in_data.p' for shift in shifts]
    nf_objs = [folder + 'Shift_' + shift + '/nf_data.p' for shift in shifts]

    nf_scrambling = scramblingGain(in_objs, nf_objs, input_method='edge', output_method='edge')

    nf_input_x = nf_scrambling[0]
    nf_input_y = nf_scrambling[1]
    nf_output_x = nf_scrambling[2]
    nf_output_y = nf_scrambling[3]
    nf_scrambling_gain = nf_scrambling[4]
    nf_input_dist = nf_scrambling[5]

    plt.figure(1)
    plt.subplot(221)
    plt.scatter(nf_input_x, nf_output_x)
    plt.xlabel('Input X [Fiber Diameter]')
    plt.ylabel('Output X [Fiber Diameter]')
    plt.subplot(222)
    plt.scatter(nf_input_y, nf_output_x)
    plt.xlabel('Input Y [Fiber Diameter]')
    plt.ylabel('Output X [Fiber Diameter]')
    plt.subplot(223)
    plt.scatter(nf_input_x, nf_output_y)
    plt.xlabel('Input X [Fiber Diameter]')
    plt.ylabel('Output Y [Fiber Diameter]')
    plt.subplot(224)
    plt.scatter(nf_input_y, nf_output_y)
    plt.xlabel('Input Y [Fiber Diameter]')
    plt.ylabel('Output Y [Fiber Diameter]')
    plt.suptitle('NF Centroid Shift')
    plt.savefig(folder + 'Near Field Shift.png')

    plt.figure(2)
    plt.title('NF Scrambling Gains')
    plt.scatter(nf_input_dist, nf_scrambling_gain)
    plt.xlabel('Input Delta [Fiber Diameter]')
    plt.ylabel('Scrambling Gain')
    plt.savefig(folder + 'Near Field SG.png')

    plt.show()

    ff_objs = [folder + 'Shift_' + shift + '/ff_data.p' for shift in shifts]

    ff_scrambling = scramblingGain(in_objs, ff_objs, input_method='edge', output_method='gaussian')

    ff_input_x = ff_scrambling[0]
    ff_input_y = ff_scrambling[1]
    ff_output_x = ff_scrambling[2]
    ff_output_y = ff_scrambling[3]
    ff_scrambling_gain = ff_scrambling[4]
    ff_input_dist = ff_scrambling[5]

    plt.figure(3)
    plt.subplot(221)
    plt.scatter(ff_input_x, ff_output_x)
    plt.xlabel('Input X [Fiber Diameter]')
    plt.ylabel('Output X [Fiber Diameter]')
    plt.subplot(222)
    plt.scatter(ff_input_y, ff_output_x)
    plt.xlabel('Input Y [Fiber Diameter]')
    plt.ylabel('Output X [Fiber Diameter]')
    plt.subplot(223)
    plt.scatter(ff_input_x, ff_output_y)
    plt.xlabel('Input X [Fiber Diameter]')
    plt.ylabel('Output Y [Fiber Diameter]')
    plt.subplot(224)
    plt.scatter(ff_input_y, ff_output_y)
    plt.xlabel('Input Y [Fiber Diameter]')
    plt.ylabel('Output Y [Fiber Diameter]')
    plt.suptitle('FF Centroid Shift')
    plt.savefig(folder + 'Far Field Shift.png')

    plt.figure(4)
    plt.title('FF Scrambling Gains')
    plt.scatter(ff_input_dist, ff_scrambling_gain)
    plt.xlabel('Input Delta [Fiber Diameter]')
    plt.ylabel('Scrambling Gain')
    plt.savefig(folder + 'Far Field SG.png')

    plt.show()