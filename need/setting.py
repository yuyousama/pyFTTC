#!/usr/bin/python3
# -*- coding:utf-8 -*-
# author:QingDoors
# 2021/8/5 21:22
import os
class Setting():
    def __init__(self):
        self.ori_img_path = ''
        self.moved_img_path = ''
        self.apply_mask = False
        self.mask_crop = False
        self.mask_img_path = ''
        self.need_algin = 1
        self.need_norm = 1
        self.correlation_method = 'circular'
        self.windowsizes = None
        self.overlap = None
        self.sig2noise_threshold = None
        self.std_factor = None
        self.median_threshold = None
        self.median_eps = None
        self.median_size = None
        self.enable_smoothn = True
        self.smoothn_p = None
        self.filter_method = None
        self.pixelsize = None
        self.young = None
        self.sigma = None
        self.lamda = None
        self.save_folder = ''
        self.save_PIV = True
        self.PIV_height = None
        self.PIV_width = None
        self.PIV_min_scale = -1
        self.PIV_max_scale = -1
        self.PIV_scale_for_lenth = 1
        self.PIV_scale_of_width = 1
        self.PIV_scale_for_arrow = 1
        self.save_FTTC = True
        self.FTTC_height = None
        self.FTTC_width = None
        self.FTTC_min_scale = -1
        self.FTTC_max_scale = -1
        self.FTTC_scale_for_lenth = 1
        self.FTTC_scale_of_width = 1
        self.FTTC_scale_for_arrow = 1
    def check(self):
        error_str = ''
        try:
            if not os.path.exists(self.ori_img_path):
                error_str += 'Image with the cell(s) doesn\'t exist.\n'
        except:
            error_str += 'Image with the cell(s) WRONG.\n'

        try:
            if not os.path.exists(self.moved_img_path):
                error_str += 'Data images don\'t exist.\n'
        except:
            error_str += 'Data images WRONG.\n'

        try:
            if self.apply_mask and (not os.path.exists(self.mask_img_path)):
                error_str += 'Mask image doesn\'t exist.\n'
        except:
            error_str += 'Mask image WRONG.\n'

        try:
            if len (self.windowsizes) != len(self.overlap):
                error_str += 'Length of Window Sizes and Overlaps should be equal.\n'
            for i in range(len(self.windowsizes)):
                if self.overlap[i] >= self.windowsizes[i]:
                    error_str += 'Overlap has to be smaller than the window_size.\n'
                    break
        except:
            error_str += 'Window Sizes and Overlaps WRONG.\n'

        try:
            if not os.path.exists(self.save_folder):
                error_str += 'Save folder doesn\'t exist.\n'
        except:
            error_str += 'Save folder WRONG.\n'

        return error_str