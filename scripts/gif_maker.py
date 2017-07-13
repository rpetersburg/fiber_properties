from fiber_properties import FiberImage, image_list
import imageio

def main():
    num_images = 1
    gap = num_images
    if num_images == 1:
        gap = 10
    max_image = 300
    cam = 'ff'
    folder = '../data/modal_noise/rv_error/coupled_ag_new/'
    ambient = image_list(folder + '../slow_agitation/ambient/' + cam + '_')
    dark = image_list(folder + '../slow_agitation/dark/' + cam + '_')
    # folder = '../data/modal_noise/rv_error/LED/'
    # folder = '../data/modal_noise/rv_error/slow_agitation/'
    # folder = '../data/modal_noise/rv_error/coupled_agitation/'
    # ambient = image_list(folder + 'ambient/' + cam + '_')
    # dark = image_list(folder + 'dark/' + cam + '_')

    with imageio.get_writer(folder + cam + '_' + str(num_images) + '.gif', mode='I') as writer:
        for i in xrange(0, max_image, gap):
            images = image_list(folder + cam + '_', num=num_images, start=i)

            image = FiberImage(images, ambient=ambient, dark=dark, camera=cam)

            print images[0]
            writer.append_data(image.get_image())

if __name__ == '__main__':
    main()
    