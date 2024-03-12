#!/usr/bin/python3
# -*- coding:utf-8 -*-
# author:QingDoors
# 2021/8/2 10:18
import numpy as np
from scipy.fft import fft2, ifft2

def FTTC(x, y, u, v, pixelsize, young, sigma, lamda):
    pixel_distance = x[0, 1] - x[0, 0]
    pixels_of_x = u.shape[1]
    pixels_of_y = v.shape[0]
    Kx_list = (creat_Ki_list(pixels_of_x) * np.pi * 2) / (pixels_of_x * pixel_distance * pixelsize)
    Ky_list = (creat_Ki_list(pixels_of_y) * np.pi * 2) / (pixels_of_y * pixel_distance * pixelsize)
    physic_u = u * pixelsize
    physic_v = v * pixelsize


    traction_force_u, traction_force_v, backup = deformation2force(physic_u, physic_v, young, sigma, Kx_list,
                                                                       Ky_list, lamda)
    return traction_force_u.real, traction_force_v.real


def deformation2force(physic_u, physic_v, young, sigma, Kx_list, Ky_list, lamda):
    '''

    :param physic_u:2d array, deformation field, u component, in uM
    :param physic_v:2d array, deformation field, v component, in uM
    :param young: float, youngs modulus, in Pa
    :param sigma: float, posson ratio
    :param Kx_list: float 1-D array, u component of wave vectors
    :param Ky_list: float 1-D array, u component of wave vectors

    :return: traction_force_u, traction_force_v, 2d array, same shape as physic_u and physic_v
    '''
    pixels_of_x = physic_u.shape[1]
    pixels_of_y = physic_u.shape[0]

    u_ft = fft2(physic_u)
    v_ft = fft2(physic_v)
    backup = (u_ft[0,0], v_ft[0,0])
    u_ft[0,0] = 0
    v_ft[0,0] = 0


    traction_force_u_ft = np.zeros_like(u_ft, dtype='complex128')
    traction_force_v_ft = np.zeros_like(v_ft, dtype='complex128')


    for j in range(pixels_of_y):
        for i in range(pixels_of_x):
            kx = Kx_list[i]
            ky = Ky_list[j]
            k = (kx ** 2 + ky ** 2) ** 0.5
            if k == 0:
                k = 1e-90

            LT = (1 - sigma) * k * k + sigma * ky * ky
            RB = (1 - sigma) * k * k + sigma * kx * kx
            if (i == (pixels_of_x / 2) + 1) or (j == (pixels_of_y / 2) + 1):
                off_diagonal = 0
            else:
                off_diagonal = -sigma * kx * ky
            A = (2 * (1 + sigma)) / (young * (k**3))
            G = np.array([[LT, off_diagonal], [off_diagonal, RB]]) * A

            GT = np.transpose(G)
            G1 = np.matmul(GT, G)
            H = np.array([[1,0],[0,1]]) * (lamda**2)
            G1 = G1 + H

            Ginv = np.linalg.inv(G1)
            uv = np.array([[u_ft[j,i]],[v_ft[j,i]]])
            GtU = np.matmul(GT, uv)
            TXY = np.matmul(Ginv, GtU)
            traction_force_u_ft[j, i] = TXY[0,0]
            traction_force_v_ft[j, i] = TXY[1,0]

    traction_force_u_ft[0,0] = 0
    traction_force_v_ft[0,0] = 0

    traction_force_u = ifft2(traction_force_u_ft)
    traction_force_v = ifft2(traction_force_v_ft)

    return traction_force_u, traction_force_v, backup


def creat_Ki_list(length):
    Ki_list = np.zeros((length,))
    i = 0
    while i <= length / 2:
        Ki_list[i] = i
        i += 1
    i = int(length / 2 - 1)+1
    while i > 0:
        Ki_list[length - i] = -i
        i -= 1
    return Ki_list
