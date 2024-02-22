# a set of functions that deal with tiff files 
# Zhilei Zhao

import tifffile as tif
import numpy as np
import matplotlib.pyplot as plt
import h5py

def ZZ_read_reshape_tiff(fileName, volume_info):
    """read single tiff file, reshape into N-D numpy array
    volume_info is a dict containing info for tif volume
    a similar matlab function exists"""
    v = tif.imread(fileName)
    # swap axis to X-Y-frames
    v = np.swapaxes(v, 0, -1)
    print('fileName:'+fileName)
    print('original dim is:', v.shape)
    # reshape into desired dimension, timePoint is inferred 
    v = np.reshape(v, (volume_info['xRes'], volume_info['yRes'], volume_info['nChannel'], volume_info['zStacks'], -1), order='F')
    print('reshaped dim is (xyczt):', v.shape)
    # check data type
    if v.dtype != volume_info['dtype']:
        print('Warning: data type is not what specified! May cause clipping!')
        print('Converted to', volume_info['dtype'])
        v = v.astype(volume_info['dtype'])
    else:
        print('data type is:', v.dtype)
    return v


def ZZ_scanimage_to_uint8(volume, saturation=0.01, mask=0.01):
    """tiff from ScanImage is of type int16
    use this function to convert it to uint8
    while adjust the constrast by saturation and mask
    a similar matlab function exists"""
    # Remove the negative part
    min_pixel_value = np.amin(volume)
    shift_volume = volume - min_pixel_value
    # mast image based saturation level
    # get the quantiles
    [low_q, high_q] = np.percentile(shift_volume, [mask*100, (1-saturation)*100])
    # clip the pixel values outside the range
    clip_volume = shift_volume
    # convert pixels greater than saturation
    clip_volume[shift_volume>high_q] = high_q
    # convert pixels smaller than mask
    clip_volume[shift_volume<low_q] = low_q
    # normalize the pixels
    normalize_volume = (clip_volume - low_q)/(high_q - low_q) * 255
    normalize_volume = normalize_volume.astype(np.dtype('uint8'))
    return normalize_volume


def ZZ_save_as_tiff(file_name, data, options={}):
    """warp tifffile imsave method"""
    # reshape into xy-frames
    data = np.reshape(data, (data.shape[0], data.shape[1], -1), order='F')
    # need to swap the axis 
    data = np.swapaxes(data, 0, -1)
    # imageJ don't support 64 bit float
    try:
        dtype = options['dtype']
    except:
        tif.imsave(file_name, data.astype(np.dtype('single')), **options)
    else:
        tif.imsave(file_name, data.astype(dtype), **options)
    print('Image save successfully as:', file_name)
    

def ZZ_running_average(movie, factor=5):
    s = movie.shape
    movie_ave = np.zeros([s[0],s[1], s[2], s[3]-factor+1], dtype='float32')
    for t in range(factor, s[3]):
        movie_ave[:,:,:,t-factor] = np.mean(movie[:,:,:,(t-factor):t], axis=3)
    return movie_ave


# define a function to pad zero stacks to a volume
def ZZ_pad_zero_volume(vol, npad=3):
    """given a 3D volume, pad zero stacks to the z axis
    vol: 3D array, size is [dim_x, dim_y, dim_z]
    npad: Int
    return: padded volumes"""
    v_size = vol.shape
    res = np.zeros((v_size[0], v_size[1], v_size[2]+2*npad), dtype=vol.dtype)
    res[:,:,npad:(npad+v_size[2])] = vol
    return(res)