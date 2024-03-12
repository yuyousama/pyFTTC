#!/usr/bin/python3
# -*- coding:utf-8 -*-
# author:QingDoors
# 2021/8/7 12:25
import numpy as np

from need.tools import crop_for_shift, align, save_data
from need.draw_from_data import draw
from need.PIV import iter_PIV
from need.FTTC import FTTC
import cv2
import os

def single_calculate(Set, img_moved, img_origin, cell_mask, file_name):

    #Step 0
    #Aligin two images and cell_mask.
    img_moved, img_origin, shift_coor = align(img_moved, img_origin)

    if Set.apply_mask:
        shift_y = shift_coor[0]
        shift_x = shift_coor[1]
        cell_mask = crop_for_shift(cell_mask, (shift_y, shift_x))
    else:
        cell_mask = 'None'

    #Step 1
    #PIV for deformation field
    x, y, u, v, mask = iter_PIV(img_moved, img_origin, Set.correlation_method, Set.need_norm, Set.windowsizes, Set.overlap,
             Set.sig2noise_threshold, Set.std_factor, Set.median_threshold, Set.median_eps, Set.median_size, Set.enable_smoothn, Set.smoothn_p, Set.filter_method)

    # Step Not_Exit
    # Crop Image With Mask

    if Set.apply_mask and Set.mask_crop:
        def get_mask(x, y, mask_img):
            mask = np.zeros_like(x, dtype=bool)
            for j in range(mask.shape[0]):
                for i in range(mask.shape[1]):
                    a = mask_img[round(y[j, i]), round(x[j, i])] != 0
                    mask[j, i] = (a)
            return mask
        grid_cell_mask = get_mask(x, y, cell_mask)
        u[grid_cell_mask == False] = 0
        v[grid_cell_mask == False] = 0

    #Step 2
    #FTTC for force field
    traction_force_u, traction_force_v = FTTC(x, y, u, v, Set.pixelsize, Set.young, Set.sigma, Set.lamda)

    #Step 3 save txt and pic
    XX = x.flatten()
    YY = y.flatten()
    UU = u.flatten()
    VV = v.flatten()
    fUU = traction_force_u.flatten()
    fVV = traction_force_v.flatten()
    save_data([XX, YY, UU, VV, fUU, fVV], os.path.join(Set.save_folder, file_name+'.txt'))

    if Set.apply_mask:
        contours, hierarchy = cv2.findContours(cell_mask, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)

    if Set.save_PIV:
        origin_pic_size = img_moved.shape
        target_pic_size = (Set.PIV_height, Set.PIV_width)
        pic, bar = draw(XX, YY, UU, VV, origin_pic_size, target_pic_size, Set.PIV_min_scale, Set.PIV_max_scale,
                          Set.PIV_scale_for_lenth, Set.PIV_scale_of_width, Set.PIV_scale_for_arrow)
        if Set.apply_mask:
            cv2.drawContours(pic, contours, -1, (255, 255, 255), 1)
        cv2.imwrite(os.path.join(Set.save_folder, file_name+'_PIV.png'), pic)
        cv2.imwrite(os.path.join(Set.save_folder, file_name+'_PIV_bar.png'), bar)

    if Set.save_FTTC:
        origin_pic_size = img_moved.shape
        target_pic_size = (Set.FTTC_height, Set.FTTC_width)
        pic, bar = draw(XX, YY, fUU, fVV, origin_pic_size, target_pic_size, Set.FTTC_min_scale, Set.FTTC_max_scale,
                          Set.FTTC_scale_for_lenth, Set.FTTC_scale_of_width, Set.FTTC_scale_for_arrow)
        if Set.apply_mask:
            cv2.drawContours(pic, contours, -1, (255, 255, 255), 1)
        cv2.imwrite(os.path.join(Set.save_folder, file_name+'_FTTC.png'), pic)
        cv2.imwrite(os.path.join(Set.save_folder, file_name+'_FTTC_bar.png'), bar)
