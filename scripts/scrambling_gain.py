from FiberProperties import ImageAnalysis, Calibration, scramblingGain
import matplotlib.pyplot as plt

folder = '../data/scrambling/2016-08-05 Prototype Core Extension 1/'

in_calibration = Calibration(dark=[folder + 'Dark/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                             ambient=[folder + 'Ambient/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

nf_calibration = Calibration(dark=[folder + 'Dark/nf_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                             ambient=[folder + 'Ambient/nf_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

ff_calibration = Calibration(dark=[folder + 'Dark/ff_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                             ambient=[folder + 'Ambient/ff_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

in_objs = []
nf_objs = []
ff_objs = []
for shift in ['00', '05', '10', '15', '20', '25', '30']:
    print 'Initializing Shift ' + shift
    in_images = [folder + 'Shift_' + shift + '/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
    nf_images = [folder + 'Shift_' + shift + '/nf_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
    ff_images = [folder + 'Shift_' + shift + '/ff_' + str(i).zfill(3) + '.fit' for i in xrange(10)]

    in_objs.append(ImageAnalysis(in_images, calibration=in_calibration, camera='in'))
    nf_objs.append(ImageAnalysis(nf_images, calibration=nf_calibration, camera='nf'))
    ff_objs.append(ImageAnalysis(ff_images, calibration=ff_calibration, camera='ff'))

nf_scrambling = scramblingGain(in_objs, nf_objs, input_method='edge', output_method='edge')

plt.figure(1)
plt.subplot(211)
plt.title('NF Centroid Shift')
plt.scatter(nf_scrambling[0], nf_scrambling[2])
plt.xlabel('Input X [um]')
plt.ylabel('Output X [um]')
plt.subplot(212)
plt.scatter(nf_scrambling[1], nf_scrambling[3])
plt.xlabel('Input Y [um]')
plt.ylabel('Output Y [um]')
plt.savefig(folder + 'Near Field Shift.png')

plt.figure(2)
plt.title('NF Scrambling Gains')
plt.plot(nf_scrambling[4])
plt.xlabel('Trial')
plt.ylabel('Scrambling Gain')
plt.savefig(folder + 'Near Field SG.png')

ff_scrambling = scramblingGain(in_objs, ff_objs, input_method='edge', output_method='gaussian')

plt.figure(3)
plt.subplot(211)
plt.title('FF Centroid Shift')
plt.scatter(ff_scrambling[0], ff_scrambling[2])
plt.xlabel('Input X [um]')
plt.ylabel('Output X [um]')
plt.subplot(212)
plt.scatter(ff_scrambling[1], ff_scrambling[3])
plt.xlabel('Input Y [um]')
plt.ylabel('Output Y [um]')
plt.savefig(folder + 'Far Field Shift.png')

plt.figure(4)
plt.title('FF Scrambling Gains')
plt.plot(ff_scrambling[4])
plt.xlabel('Trial')
plt.ylabel('Scrambling Gain')
plt.savefig(folder + 'Far Field SG.png')

plt.show()