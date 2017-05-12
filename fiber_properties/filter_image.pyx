"""filter_image.pyx was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The functions in this module attempt to speed up median filtering using cython
"""
from __future__ import division
import numpy as np
cimport numpy as np
import math
from scipy.signal import medfilt2d

def median(np.ndarray[np.float_t] array):
    cdef int length = len(array)
    cdef int pivot
    cdef object kth

    if length % 2 == 0:
        pivot = length // 2
        kth = [pivot-1, pivot]
        array.partition(kth)
        return (array[pivot-1] + array[pivot]) / 2.0

    pivot = (length - 1) // 2
    kth = [pivot]
    array.partition(kth)
    return array[pivot]

# def median(np.ndarray array): # only for odd length arrays
#     cdef int pivot = (len(array) - 1) // 2
#     array.partition(pivot)
#     return array[pivot]

# cdef float _median(np.ndarray array):
#     cdef int length = len(array)
#     cdef bool even = length % 2 == 0
#     cdef int pivot = length // 2
#     cdef int kth
#     if even:
#         kth = [pivot-1, pivot]
#     else:
#         kth = [pivot]
#     array.partition(kth)

#     if even:
#         return (array[pivot-1] + array[pivot]) / 2.0
#     return array[pivot]

def c_filter_image(np.ndarray image, int kernel_size):
    cdef int height, width
    height = int(image.shape[0])
    width = int(image.shape[1])

    cdef int radius = (kernel_size-1) // 2
    cdef np.ndarray mask, temp_mask, x_array, y_array

    x_array, y_array = np.meshgrid(np.arange(kernel_size),
                                   np.arange(kernel_size))
    mask = (x_array-radius)**2 + (y_array-radius)**2 <= radius**2

    cdef np.ndarray image_crop, filtered_image
    filtered_image = np.zeros_like(image)

    cdef int x, y, top, bottom, left, right, mask_top, mask_bottom, mask_left, mask_right
    for y in xrange(height):
        if y < radius:
            mask_top = kernel_size - radius - 1 - y
            top = 0
        else:
            mask_top = 0
            top = y - radius
        if y > height - 1 - radius:
            mask_bottom = kernel_size - radius + (height - 1 - y)
            bottom = height
        else:
            mask_bottom = kernel_size
            bottom = y + radius + 1

        for x in xrange(width):
            if x < radius:
                mask_left = kernel_size - radius - 1 - x
                left = 0
            else:
                mask_left = 0
                left = x - radius
            if x > width - 1 - radius:
                mask_right = kernel_size - radius + (width - 1 - x)
                right = width
            else:
                mask_right = kernel_size
                right = x + radius + 1

            temp_mask = mask[mask_top:mask_bottom, mask_left:mask_right]
            image_crop = image[top:bottom, left:right]
            filtered_image[y, x] = median(image_crop[temp_mask])

    return filtered_image

def c_filter_image_zero_fill(np.ndarray image, int kernel_size):
    cdef int height, width
    height = int(image.shape[0])
    width = int(image.shape[1])

    cdef int radius = (kernel_size-1) // 2
    cdef np.ndarray mask, x_array, y_array

    x_array, y_array = np.meshgrid(np.arange(kernel_size),
                                   np.arange(kernel_size))
    mask = (x_array-radius)**2 + (y_array-radius)**2 <= radius**2

    cdef np.ndarray image_crop, filtered_image, zero_image
    filtered_image = np.zeros_like(image)
    zero_image = np.zeros([kernel_size, kernel_size])
    cdef int x, y, top, bottom, left, right, mask_top, mask_bottom, mask_left, mask_right
    for y in xrange(height):
        if y < radius:
            mask_top = kernel_size - radius - 1 - y
            top = 0
        else:
            mask_top = 0
            top = y - radius
        if y > height - 1 - radius:
            mask_bottom = kernel_size - radius + (height - 1 - y)
            bottom = height
        else:
            mask_bottom = kernel_size
            bottom = y + radius + 1

        for x in xrange(width):
            if x < radius:
                mask_left = kernel_size - radius - 1 - x
                left = 0
            else:
                mask_left = 0
                left = x - radius
            if x > width - 1 - radius:
                mask_right = kernel_size - radius + (width - 1 - x)
                right = width
            else:
                mask_right = kernel_size
                right = x + radius + 1

            image_crop = zero_image.copy()
            image_crop[mask_top:mask_bottom, mask_left:mask_right] = image[top:bottom, left:right]
            filtered_image[y, x] = median(image_crop[mask])

    return filtered_image

# cdef float[:,:] _subset(float[:,:] array, int top, int bottom, int left, int right):
#     cdef float[bottom-top, right-left] subframe
#     cdef int i, j
#     for i in xrange(right-left):
#         for j in xrange(bottom-top):
#             subframe[j][i] = array[top+j][left+i]
#     return subframe

# def _median(np.ndarray array):
#     return _select(array, int((array.size-1)/2))

# def _select(np.ndarray array, int k):
#     cdef int array_size = array.size
#     if array_size < 6:
#         return float(np.sort(array)[k])
#     cdef np.ndarray[np.float_t] medians = np.empty(int(math.ceil(array_size/5.0)), dtype=np.float)
#     cdef int i
#     for i in xrange(0, array_size, 5):
#         medians[int(i/5)] = _median(np.sort(array[i:i+5]))
#     cdef float median = _median(medians)

#     cdef np.ndarray lower = array[array<median]
#     cdef np.ndarray upper = array[array>median]
#     cdef int lower_size = lower.size
#     if k < lower_size:
#         return _select(lower, k)
#     elif k > lower_size:
#         return _select(upper, k-lower_size-1)
#     else:
#         return median

# def median(np.ndarray[float, mode='c'] array):
#     # cdef float[:] y = array
#     # cdef float x[]
#     # x = y
#     # return _median(array.shape[0], &array[0])



# cdef float _median(int n, float *x):
#     cdef float temp
#     cdef int i, j
#     for i in xrange(n):
#         for j in xrange(i, n):
#             if x[j] < x[i]:
#                 temp = x[i]
#                 x[i] = x[j]
#                 x[j] = temp

#     if n % 2 == 0:
#         return (x[n/2] + x[n/2 - 1]) / 2.0
#     return x[n/2]