#!/usr/bin/python3
# -*- coding:utf-8 -*-
# author:QingDoors
# 2021/8/1 16:08
import numpy as np
from openpiv.pyprocess import extended_search_area_piv, get_coordinates, get_field_shape
from openpiv import smoothn
from scipy.interpolate import RectBivariateSpline
import scipy.ndimage as scn

from need.tools import typical_validation
def iter_PIV(image_a, image_b, correlation_method, normalized_correlation, windowsizes, overlap,
             sig2noise_threshold, std_factor, median_threshold, median_eps, median_size, enable_smoothn, smoothn_p, filter_method):
    '''

    :param image_a: 2d array
    :param image_b: 2d array
    :param correlation_method: string, 'circular' or 'linear', latter is slower and requires also normalized_correlation = True
    :param normalized_correlation: bool, if True, then the image intensity will be modified by removing
        the mean, dividing by the standard deviation and the correlation map will be normalized. It's slower but could be more robust.
    :param windowsizes: tuple, window_size of each pass, everyone should be a power of 2
    :param overlap: turple, overlap of each pass, same length of windowsizes. 50% of windowsizes recommend.
    :param sig2noise_threshold: None or float. 1.05 recommend. If none, sig2noise validation will be disabled.
    :param std_factor: None or float. 10.0 recommend. If none, STD validation will be disabled.
    :param median_threshold: None or float. 2.0 recommend. If none, NMT validation will be disabled.
    :param median_eps: float. available only when median_threshold is float. 0.1 recommend
    :param median_size: int. available only when median_threshold is float. 1 recommend
    :param enable_smoothn: bool, Enables smoothing.
    :param smoothn_p: None or float, parameter of smoothn. 0.05 recommend, or set it to None so it will be calculated automatically
    :param filter_method: string, method to replace outliers. 'localmean', 'disk' or 'distance'
    :return:x,y,u,v,masK: 2d np.ndarray
    '''
    x, y, u, v, s2n = first_pass(
        image_a,
        image_b,
        window_size=windowsizes[0],
        overlap=overlap[0],
        correlation_method = correlation_method,
        normalized_correlation = normalized_correlation,
    )

    u, v, mask = typical_validation(u, v, s2n, sig2noise_threshold, std_factor, median_threshold, median_eps, median_size, filter_method)
    u = np.ma.masked_array(u, mask=np.ma.nomask)
    v = np.ma.masked_array(v, mask=np.ma.nomask)
    if enable_smoothn:
        u, dummy_u1, dummy_u2, dummy_u3 = smoothn.smoothn(u, s=smoothn_p)
        v, dummy_v1, dummy_v2, dummy_v3 = smoothn.smoothn(v, s=smoothn_p)

    u = np.ma.masked_array(u, mask=np.ma.nomask)
    v = np.ma.masked_array(v, mask=np.ma.nomask)

    for i in range(1, len(windowsizes)):
        x, y, u, v, s2n, mask = multipass_img_deform(
        image_a,
        image_b,
        x,
        y,
        u,
        v,
        windowsizes[i],
        overlap[i],
        correlation_method,
        normalized_correlation,
        sig2noise_threshold,
        std_factor,
        median_threshold,
        median_eps,
        median_size,
        filter_method)
        if enable_smoothn and i < len(windowsizes) - 1:
            u, dummy_u1, dummy_u2, dummy_u3 = smoothn.smoothn(u, s=smoothn_p)
            v, dummy_v1, dummy_v2, dummy_v3 = smoothn.smoothn(v, s=smoothn_p)
        u = np.ma.masked_array(u, np.ma.nomask)
        v = np.ma.masked_array(v, np.ma.nomask)
    u = u.filled(0.)
    v = v.filled(0.)
    u = np.ma.masked_array(u, np.ma.nomask)
    v = np.ma.masked_array(v, np.ma.nomask)

    return x,y,u,v,mask



def first_pass(image_a, image_b, window_size, overlap, correlation_method, normalized_correlation):
    u, v, s2n = extended_search_area_piv(
        image_a,
        image_b,
        window_size=window_size,
        overlap=overlap,
        search_area_size=window_size,
        width=2,
        subpixel_method='gaussian',
        sig2noise_method='peak2peak',
        correlation_method=correlation_method,
        normalized_correlation=normalized_correlation
    )

    shapes = np.array(get_field_shape(image_a.shape, window_size, overlap))
    u = u.reshape(shapes)
    v = v.reshape(shapes)
    s2n = s2n.reshape(shapes)

    x, y = get_coordinates(image_a.shape, window_size, overlap)

    return x, y, u, v, s2n

def multipass_img_deform(image_a, image_b, x_old, y_old, u_old, v_old, window_size, overlap, correlation_method, normalized_correlation,
        sig2noise_threshold, std_factor, median_threshold, median_eps, median_size, filter_method):
    x, y = get_coordinates(image_a.shape, window_size, overlap)

    y_old = y_old[:, 0]
    x_old = x_old[0, :]

    y_int = y[:, 0]
    x_int = x[0, :]

    ip = RectBivariateSpline(y_old, x_old, u_old.filled(0.))
    u_pre = ip(y_int, x_int)

    ip2 = RectBivariateSpline(y_old, x_old, v_old.filled(0.))
    v_pre = ip2(y_int, x_int)

    x_new, y_new, ut, vt = create_deformation_field(image_a, x, y, u_pre, v_pre)
    image_a = scn.map_coordinates(image_a, ((y_new - vt / 2, x_new - ut / 2)),order=5, mode='nearest')
    image_b = scn.map_coordinates(image_b, ((y_new + vt / 2, x_new + ut / 2)),order=5, mode='nearest')


    u, v, s2n = extended_search_area_piv(
        image_a,
        image_b,
        window_size=window_size,
        overlap=overlap,
        search_area_size=window_size,
        width=2,
        subpixel_method="gaussian",
        sig2noise_method="peak2peak",
        correlation_method=correlation_method,
        normalized_correlation=normalized_correlation,
    )

    shapes = np.array(get_field_shape(image_a.shape,
                                      window_size,
                                      overlap))
    u = u.reshape(shapes)
    v = v.reshape(shapes)
    s2n = s2n.reshape(shapes)

    u += u_pre
    v += v_pre

    u, v, mask = typical_validation(u, v, s2n, sig2noise_threshold, std_factor, median_threshold, median_eps, median_size, filter_method)

    u = np.ma.masked_array(u, np.ma.nomask)
    v = np.ma.masked_array(v, np.ma.nomask)

    return x, y, u, v, s2n, mask

def create_deformation_field(frame, x, y, u, v, kx=3, ky=3):
    y1 = y[:, 0]
    x1 = x[0, :]
    side_x = np.arange(frame.shape[1])
    side_y = np.arange(frame.shape[0])
    ip = RectBivariateSpline(y1, x1, u, kx=kx, ky=ky)
    ut = ip(side_y, side_x)
    ip2 = RectBivariateSpline(y1, x1, v, kx=kx, ky=ky)
    vt = ip2(side_y, side_x)

    x, y = np.meshgrid(side_x, side_y)
    return x, y, ut, vt
