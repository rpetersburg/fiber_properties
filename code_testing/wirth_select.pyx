import pyximport; pyximport.install()
# import numpy as np
# cimport numpy as np
# ctypedef np.double_t T # or whatever

import numpy as np

def median(x):

    """ median in average O(n) time """

    n = x.shape[0]
    k = n >> 1
    s = wirthselect(x, k)
    if n & 1:
        return s[k]
    else:
        return 0.5*(s[k]+s[:k].max())

def wirthselect(array, k):
   
    """ Niklaus Wirth's selection algortithm """

    a = np.ascontiguousarray(array)
    if (a is array): a = a.copy()

    l = 0
    m = a.shape[0] - 1
    while l < m:
        x = a[k]
        i = l
        j = m
        while 1:
            while a[i] < x: i += 1
            while x < a[j]: j -= 1
            if i <= j:
                tmp = a[i]
                a[i] = a[j]
                a[j] = tmp
                i += 1
                j -= 1
            if i > j: break
        if j < k: l = i
        if k < i: m = j

    return a

# def wirthselect(np.ndarray[T, ndim=1, mode='c'] array, int k):
   
#     cdef int i, j, l, m
#     cdef T x, tmp
#     cdef T* a

#     _array = np.ascontiguousarray(array)
#     if (_array is array): _array = _array.copy()
#     a = <T*> _array.data

#     l = 0
#     m = <int> a.shape[0] - 1
#     with nogil:
#         while l < m:
#             x = a[k]
#             i = l
#             j = m
#             while 1:
#                 while a[i] < x: i += 1
#                 while x < a[j]: j -= 1
#                 if i <= j:
#                     tmp = a[i]
#                     a[i] = a[j]
#                     a[j] = tmp
#                     i += 1
#                     j -= 1
#                 if i > j: break
#             if j < k: l = i
#             if k < i: m = j

#     return _array