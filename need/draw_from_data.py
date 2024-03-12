#!/usr/bin/python3
# -*- coding:utf-8 -*-
# author:QingDoors
# 2021/8/5 16:06
import numpy as np
import cv2
import os




def draw_triangle(img, p, l, vector, c):
    '''
    :param img:图片，array
    :param p: 底边中点坐标,int
    :param l:三角形高,float/int
    :param vector:高的向量,float
    :param c:颜色
    :return:img
    '''
    u = vector[0]
    v = vector[1]
    lengyh_of_vector = (u**2 + v**2)**0.5


    p0 = p
    x0 = p0[0]
    y0 = p0[1]

    x1 =int(round(x0 + (l * u) / lengyh_of_vector))
    y1 =int(round(y0 + (l * v) / lengyh_of_vector))
    p1 = (x1, y1)

    d = l * np.tan((15/180) * np.pi)
    v_vector = (v / lengyh_of_vector, -u/lengyh_of_vector)

    v_u = d * v_vector[0]
    v_v = d * v_vector[1]

    x2 = int(round(x0 + v_u))
    y2 = int(round(y0 + v_v))
    p2 = (x2, y2)

    x3 = int(round(x0 - v_u))
    y3 = int(round(y0 - v_v))
    p3 = (x3, y3)

    pts = np.array([p1, p2, p3], np.int32)
    cv2.fillConvexPoly(img, pts, c)  # 可以填充
    return (img)
def draw_arrow(img, p, vector, c, width, min_arrow):
    '''

    :param img: 背景图片，array
    :param p: 起始点坐标，（x,y）,int
    :param vector: 向量，(u, v),float
    :param c: 颜色
    width:线宽,int
    min_arrow,最小箭头长度，int
    :return:img(draw)
    '''
    p1 = p
    p2 = (int(round(p[0]+vector[0])), int(round(p[1]+vector[1])))
    cv2.line(img, p1, p2, c, width)
    length = (vector[0] ** 2 + vector[1] ** 2) ** 0.5
    l = 0.8*length
    if l < min_arrow:
        l = min_arrow
    if vector[0] != 0 or vector[1] != 0:
        img = draw_triangle(img, p2, l, vector, c)
    return (img)
def scale_data(data, f_min=None, f_max=None):
    k = 499.0 / (f_max - f_min)
    b = (499.0 * f_min) / (f_min - f_max)

    result = data * k
    result = result + b
    return (result, k, b)

def convert(number, decimal=2):
    if number * np.power(10, decimal) < 1:
        e = 0
        a = number
        while True:
            e += 1
            a *= 10
            if a>=1:
                break
        return str(round(a, decimal)) + 'e-' + str(e)
    else:
        return str(round(number, decimal))


eg_bar = cv2.imread(f'{os.getcwd()}/need/eg_bar.png')
def creat_bar(min_scale, max_scale, size):
    width = max(round(size/15.0), 20)
    bar = np.zeros((width +50, size, 3))
    for i in range(size):
        c = eg_bar[round((i*499)/size), 0, :].astype('float64')
        bar[:width,size-i-1,0] = c[0]
        bar[:width,size-i-1,1] = c[1]
        bar[:width,size-i-1,2] = c[2]

    bar = cv2.putText(bar, convert(max_scale, 2), (size-125, width+40), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255,255,255))
    bar = cv2.putText(bar, convert(min_scale, 2), (0, width+40), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255,255,255))
    return bar


def draw(XX, YY, UU, VV, origin_pic_size, target_pic_size, min_scale_bar, max_scale_bar, scale_for_lenth=1, scale_of_width=1, scale_for_arrow=1):
    distance = round(XX[1] - XX[0])
    if target_pic_size == (None, None):
        target_pic_size = origin_pic_size

    ZZ = np.power(UU, 2) + np.power(VV, 2)
    ZZ = np.sqrt(ZZ)
    if min_scale_bar < 0:
        min_force = np.min(ZZ)
        max_force = np.max(ZZ)
    else:
        min_force = min_scale_bar
        max_force = max_scale_bar

    print(format(min_scale_bar, '.2f'), format(max_scale_bar, '.2f'))
    print(format(min_force, '.2f'), format(max_force, '.2f'))

    ZZ, k, b = scale_data(ZZ, min_force, max_force)

    pic = np.zeros((target_pic_size[0], target_pic_size[1], 3))
    
    UU = UU * k + b
    VV = VV * k + b
    ZZ = ZZ * -1 + 499

    ZZ = np.round(ZZ)
    ZZ[ZZ > 499] = 499
    ZZ[ZZ < 0] = 0

    UU = UU * scale_for_lenth
    VV = VV * scale_for_lenth

    scale_for_x = target_pic_size[1] / origin_pic_size[1]
    scale_for_y = target_pic_size[0] / origin_pic_size[0]
    target_distance = distance * min(scale_for_y, scale_for_x)
    width = max(round(target_distance * 0.1 * scale_of_width), 1)
    min_arrow = max(round(target_distance * 0.3 * scale_for_arrow), 1)
    for i in range(len(XX)):
        try:
            x = round(XX[i] * scale_for_x)
            y = round(YY[i] * scale_for_y)
            u = UU[i]
            v = VV[i]
            z = int(ZZ[i])
            u = (u * distance * scale_for_x) / (499 * 1.8)
            v = (v * distance * scale_for_y) / (499 * 1.8)
            c = eg_bar[z, 0, :].astype('float64')
            pic = draw_arrow(pic, (x, y), (u, v), c, width, min_arrow)
        except:
            pass
    bar = creat_bar(min_force, max_force, max(target_pic_size))
    return pic, bar






