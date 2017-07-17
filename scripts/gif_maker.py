from fiber_properties import FiberImage, image_list
import imageio

def main():
    num_images = 1
    interval = 1
    max_image = 10
    cam = 'nf'
    folder = '../data/modal_noise/Kris_data/octagonal_200um/coupled_agitation/'
    # folder = '../data/modal_noise/rv_error/coupled_agitation/'
    ambient = image_list(folder + '../ambient/' + cam + '_')
    dark = image_list(folder + '../dark/' + cam + '_')

    # with imageio.get_writer(folder + cam + '_' + str(num_images) + '.gif', mode='I') as writer:
    with imageio.get_writer(folder + cam + '_snr.gif', mode='I') as writer:
        for i in xrange(0, max_image, interval):
            images = image_list(folder + cam + '_', num=i+1, start=0)

            image = FiberImage(images, ambient=ambient, dark=dark, camera=cam)

            print images[-1]
            writer.append_data(image.get_image())

if __name__ == '__main__':
    main()
    