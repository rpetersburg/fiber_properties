from FiberProperties import scramblingGain, imageList
import matplotlib.pyplot as plt

if __name__ == '__main__':

    NEW_DATA = False
    folder = '../data/scrambling/2016-08-05 Prototype Core Extension 1/'

    if NEW_DATA:
        from FiberProperties import ImageAnalysis

        in_dark = imageList(folder + 'Dark/in_')
        in_ambient = imageList(folder + 'Ambient/in_')
        nf_dark = imageList(folder + 'Dark/nf_')
        nf_ambient = imageList(folder + 'Ambient/nf_')
        ff_dark = imageList(folder + 'Dark/ff_')
        ff_ambient = imageList(folder + 'Ambient/ff_')

        for shift in ['00', '05', '10', '15', '20', '25', '30']:
            print 'Initializing Shift ' + shift
            in_images = imageList(folder + 'Shift_' + shift + '/in_')
            nf_images = imageList(folder + 'Shift_' + shift + '/nf_')
            ff_images = imageList(folder + 'Shift_' + shift + '/ff_')

            ImageAnalysis(in_images, in_dark, in_ambient, camera='in').save()
            ImageAnalysis(nf_images, nf_dark, nf_ambient, camera='nf').save()
            ImageAnalysis(ff_images, ff_dark, ff_ambient, camera='ff').save()

    shifts = ['00', '05', '10', '15', '20', '25', '30']
    in_objs = [folder + 'Shift_' + shift + '/in_object.pkl' for shift in shifts]
    nf_objs = [folder + 'Shift_' + shift + '/nf_object.pkl' for shift in shifts]

    nf_scrambling = scramblingGain(in_objs, nf_objs, input_method='edge', output_method='edge')

    nf_input_x = nf_scrambling[0]
    nf_input_y = nf_scrambling[1]
    nf_output_x = nf_scrambling[2]
    nf_output_y = nf_scrambling[3]
    nf_scrambling_gain = nf_scrambling[4]
    nf_input_dist = nf_scrambling[5]
    nf_output_dist = nf_scrambling[6]

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
    plt.scatter(nf_input_dist, nf_output_dist)
    plt.xlabel('Input Delta [Fiber Diameter]')
    plt.ylabel('Output Delta [Fiber Diameter]')
    plt.savefig(folder + 'Near Field SG.png')

    plt.show()

    ff_objs = [folder + 'Shift_' + shift + '/ff_object.pkl' for shift in shifts]

    ff_scrambling = scramblingGain(in_objs, ff_objs, input_method='edge', output_method='gaussian')

    ff_input_x = ff_scrambling[0]
    ff_input_y = ff_scrambling[1]
    ff_output_x = ff_scrambling[2]
    ff_output_y = ff_scrambling[3]
    ff_scrambling_gain = ff_scrambling[4]
    ff_input_dist = ff_scrambling[5]
    ff_output_dist = ff_scrambling[6]

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
    plt.scatter(ff_input_dist, ff_output_dist)
    plt.xlabel('Input Delta [Fiber Diameter]')
    plt.ylabel('Output Delta [Fiber Diameter]')
    plt.savefig(folder + 'Far Field SG.png')

    plt.show()