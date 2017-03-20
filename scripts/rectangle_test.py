from fiber_properties import ImageAnalysis, image_list

obj = ImageAnalysis(image_list('pos_1/nf_'),
                    ambient=image_list('../ambient/nf_'),
                    camera='ff')

print obj.get_fiber_center(method='rectangle', show_image=True)