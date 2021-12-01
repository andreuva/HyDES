import numpy as np

def deriv2D(array, dx=None, dy=None, axis=0):

    if axis == 0:
        deriv = ( (array[:-1,1:] - array[:-1,:-1])/dx + (array[1:,1:] - array[1:,:-1])/dx )/2.0
    elif axis == 1:
        deriv = ( (array[1:,:-1] - array[:-1,:-1])/dy + (array[1:,1:] - array[:-1,1:])/dy )/2.0
    else:
        raise ValueError("Must specify axis = 0 (x) or axis = 1 (y)")

    return deriv


def midval(array):
    return (array[:-1,:-1] + array[1:,:-1] + array[:-1,1:] + array[1:,1:])/4.0