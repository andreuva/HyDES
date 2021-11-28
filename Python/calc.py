import numpy as np

def deriv2D(array, dx=None, dy=None):

    if dx:
        deriv = ( (array[1:,0:-1] - array[0:-1,0:-1])/dx + (array[1:,1:] - array[0:-1,1:])/dx )/2.0
    elif dy:
        deriv = ( (array[0:-1,1:] - array[0:-1,0:-1])/dy + (array[1:,1:] - array[1:,0:-1])/dy )/2.0
    else:
        raise ValueError("Must specify dx or dy")

    return deriv


def midval(array):
    return (array[:-1,:-1] + array[1:,:-1] + array[:-1,1:] + array[1:,1:])/4.0