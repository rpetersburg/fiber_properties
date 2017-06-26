from fiber_properties import FiberImage, image_list
import imageio

num_images = 1
max_image = 300
cam = 'nf'
# folder = '../data/modal_noise/rv_error/coupled_agitation/'
# folder = '../data/modal_noise/rv_error/LED/'
folder = '../data/modal_noise/rv_error/slow_agitation/'
# folder = '../data/modal_noise/amp_freq_200um/unagitated_1s'

with imageio.get_writer(folder + cam + '_' + str(num_images) + '.gif', mode='I') as writer:
    for i in xrange(0, max_image, num_images):
        images = image_list(folder + cam + '_', num=num_images, start=i)
        ambient = image_list(folder + 'ambient/' + cam + '_', num=num_images)
        dark = image_list(folder + 'dark/' + cam + '_', num=num_images)

        image = FiberImage(images, ambient=ambient, dark=dark, camera=cam)

        print images[0]
        writer.append_data(image.get_image())