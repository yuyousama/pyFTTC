#!/usr/bin/python3
# -*- coding:utf-8 -*-
# author:QingDoors
# 2021/7/31 22:28

import cv2
from scipy.ndimage.interpolation import shift
from skimage.registration import phase_cross_correlation
import numpy as np
from openpiv.filters import replace_outliers

def crop_for_shift(img, shift_coor):
    shift_y = round(shift_coor[0])
    shift_x = round(shift_coor[1])
    if shift_x <= 0:
        img = img[:, -shift_x:]
    else:
        img = img[:, :-shift_x]
    if shift_y <= 0:
        img = img[-shift_y:, :]
    else:
        img = img[:-shift_y, :]
    return img

def align(image1, image2):
    '''
    Align 2 images. Only image1 will be shifted.
    :param image1: 2d array. First image to be aligned. Will be shifted and croped.
    :param image2: 2d array. Second image to be aligned. Will be croped.
    :return:
    image1:Aligned image1, shifted and croped.
    image2:Aligned image2, croped.
    shift_coor: Shift vector of the 2 images, maybe useful later.
    '''
    shift_coor = phase_cross_correlation(image1, image2, upsample_factor=100)
    shift_y = shift_coor[0][0]
    shift_x = shift_coor[0][1]
    image1 = shift(image1, shift=(-shift_y, -shift_x), order=5)
    image1 = crop_for_shift(image1, (shift_y, shift_x))
    image2 = crop_for_shift(image2, (shift_y, shift_x))
    return image1, image2, (shift_y, shift_x)

def typical_validation(u, v, s2n, sig2noise_threshold, std_factor, median_threshold, median_eps, median_size, filter_method):
    mask = np.zeros(u.shape, dtype=bool)

    if sig2noise_threshold != None:
        mask[s2n < sig2noise_threshold] = True

    if std_factor != None:
        z = u ** 2 + v ** 2
        m = np.nanmean(z)
        std_threshold = m + np.nanstd(z) * std_factor
        mask[z>std_threshold] = True
    if median_threshold != None:
        for i in np.arange(median_size, u.shape[0]-median_size, 1):
            for j in np.arange(median_size, u.shape[1]-median_size, 1):
                if mask[i,j]:
                    continue
                neighbor_u = u[i-median_size:i+median_size+1, j-median_size:j+median_size+1].flatten()
                neighbor_v = v[i-median_size:i+median_size+1, j-median_size:j+median_size+1].flatten()

                neighbor_u = np.append(neighbor_u[:2*median_size*(median_size+1)], neighbor_u[2*median_size*(median_size+1)+1:])
                neighbor_v = np.append(neighbor_v[:2*median_size*(median_size+1)], neighbor_v[2*median_size*(median_size+1)+1:])

                median_u = np.nanmedian(neighbor_u)
                median_v = np.nanmedian(neighbor_v)

                fluct_u = u[i,j] - median_u
                fluct_v = v[i,j] - median_v

                res_u = neighbor_u - median_u
                res_v = neighbor_v - median_v

                medianRes_u = np.nanmedian(abs(res_u))
                medianRes_v = np.nanmedian(abs(res_v))

                normFluct_u = abs(fluct_u / (medianRes_u + median_eps))
                normFluct_v = abs(fluct_v / (medianRes_v + median_eps))
                normFluct = (normFluct_u ** 2 + normFluct_v ** 2) ** 0.5
                if normFluct > median_threshold:
                    mask[i,j] = True
    u[mask] = np.nan
    v[mask] = np.nan
    u, v = replace_outliers(
        u,
        v,
        method=filter_method,
        max_iter=10,
        kernel_size=2
    )
    return u, v, mask
def save_data(data, file):
    out = np.vstack([m.flatten() for m in data])
    np.savetxt(
        file,
        out.T,
        delimiter='\t',
        header="x"
        + '\t'
        + "y"
        + '\t'
        + "u"
        + '\t'
        + "v"
        + '\t'
        + "f_u"
        + '\t'
        + 'f_v',
    )
