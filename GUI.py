from tkinter import ttk
import tkinter.messagebox
from tkinter.filedialog import *
import numpy as np
import cv2
import os
import threading

from need.setting import Setting
from need.single_calculate import single_calculate
Set = Setting()

def check_parameters():
    err_str = Set.check()
    if len(err_str) > 0:
        tkinter.messagebox.showerror(title='ERROR！', message=err_str)
    else:
        strat_button['state'] = 'disabled'
        __calculate()
def __calculate():
    T = threading.Thread(target=calculate, args=())
    T.setDaemon(True)
    T.start()
def update_pb(n, size):
    a = int((n*100)/size)
    pb.set(a)

def calculate():
    img_origin = cv2.imread(Set.ori_img_path)[:,:,0]
    if Set.apply_mask:
        cell_mask = cv2.imread(Set.mask_img_path)[:,:,0]
    else:
        cell_mask = None
    img_moved_list = os.listdir(Set.moved_img_path)
    error_list = []
    pb.set(0)
    n = 0
    import time

    for file_name in img_moved_list:
        a = time.time()
        img_moved = cv2.imread(os.path.join(Set.moved_img_path, file_name))[:,:,0]
        single_calculate(Set, img_moved, img_origin, cell_mask, file_name[:-4])
        print (time.time() - a)
        # try:
        #     single_calculate(Set, img_moved, img_origin, cell_mask, file_name[:-4])
        # except:
        #     error_list.append(file_name)
        n += 1
        update_pb(n, len(img_moved_list))
    if len(error_list) == 0:
        tkinter.messagebox.showinfo(title='Finish!', message='Mission Complete Successfully.')
    else:
        error_msg = 'Mission Complete but something WRONG with:'
        for name in error_list:
            error_msg += '\n'
            error_msg += name
        tkinter.messagebox.showwarning(title='Finish!', message=error_msg)
    strat_button['state'] = 'normal'

#----------------------------------------------------------------------------------------------------------#
root = Tk()
root.title("pyFTTC  By NJU Zhang Yu")
root.grid_columnconfigure(0, weight=1)
root.grid_columnconfigure(1, weight=1)
root.grid_rowconfigure(0, weight=5)
root.grid_rowconfigure(1, weight=10)
#----------------------------------------------------------------------------------------------------------#
mainframe_0 = ttk.Frame(root, padding="3 3 12 12")
mainframe_0.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe_0['borderwidth'] = 2
mainframe_0['relief'] = 'raised'

mainframe_0.grid_columnconfigure(0, weight=1)
# mainframe_0.grid_columnconfigure(1, weight=1)
mainframe_0.grid_rowconfigure(0, weight=1)
mainframe_0.grid_rowconfigure(1, weight=1)
mainframe_0.grid_rowconfigure(2, weight=1)
mainframe_0.grid_rowconfigure(3, weight=1)
mainframe_0.grid_rowconfigure(4, weight=1)


ttk.Label(mainframe_0, text="1.Choose your data files.").grid(column=0, row=0, padx=5, pady=5, sticky = W)

def set_withcell_file():
    filepath = askopenfilename()
    Set.ori_img_path = filepath
    withcell_file_name.set(filepath)

withcell_file_name = StringVar()
ttk.Entry(mainframe_0, textvariable=withcell_file_name, state='readonly').grid(row=1, column=0, padx=5, pady=5, sticky = (W,E))
withcell_file_name.set('Image with the cell(s).')

ttk.Button(mainframe_0, text = 'Open...', command = set_withcell_file).grid(row=1, column=1, padx=5, pady=5, sticky = W)


def set_without_folder():
    folderpath = askdirectory()
    Set.moved_img_path = folderpath
    withoutcell_folder_name.set(folderpath)
withoutcell_folder_name = StringVar()
ttk.Entry(mainframe_0, textvariable=withoutcell_folder_name, state='readonly' ).grid(row=2, column=0, padx=5, pady=5, sticky = (W,E))
withoutcell_folder_name.set('Folder of images.')

ttk.Button(mainframe_0, text = 'Open...', command = set_without_folder).grid(row=2, column=1, padx=5, pady=5, sticky = W)

# def set_apply_mask():
#     Set.apply_mask = apply_mask.get() == '1'
#
# apply_mask = StringVar()
# ttk.Checkbutton(mainframe_0, text='Apply Mask', command=set_apply_mask,variable=apply_mask).grid(column=0, row=3, padx=5, pady=5, sticky = W)
# apply_mask.set(0)
#
#
# def set_mask_crop():
#     Set.mask_crop = mask_crop.get() == '1'
#
#
# mask_crop = StringVar()
# ttk.Checkbutton(mainframe_0, text='Crop Outside', command=set_mask_crop,variable=mask_crop).grid(column=1, row=3, padx=5, pady=5, sticky = W)
# mask_crop.set(0)
#
#
#
# def set_mask_file():
#     filepath = askopenfilename()
#     Set.mask_img_path = filepath
#     mask_file_name.set(filepath)
#
# mask_file_name = StringVar()
# ttk.Entry(mainframe_0, textvariable=mask_file_name,  state='readonly').grid(row=4, column=0, padx=5, pady=5, sticky = (W,E))
# mask_file_name.set('Choose mask file.')
#
# ttk.Button(mainframe_0, text = 'Open...', command = set_mask_file).grid(row=4, column=1, padx=5, pady=5, sticky = W)

#----------------------------------------------------------------------------------------------------------#

mainframe_1 = ttk.Frame(root, padding="3 3 12 12")
mainframe_1.grid(column=0, row=1, sticky=(N, W, E, S))
mainframe_1['borderwidth'] = 2
mainframe_1['relief'] = 'raised'

mainframe_1.grid_columnconfigure(0, weight=1)
mainframe_1.grid_columnconfigure(1, weight=1)
mainframe_1.grid_rowconfigure(0, weight=1)
mainframe_1.grid_rowconfigure(1, weight=1)
mainframe_1.grid_rowconfigure(2, weight=1)
mainframe_1.grid_rowconfigure(3, weight=1)
mainframe_1.grid_rowconfigure(4, weight=1)
mainframe_1.grid_rowconfigure(5, weight=1)
mainframe_1.grid_rowconfigure(6, weight=1)
mainframe_1.grid_rowconfigure(7, weight=1)
mainframe_1.grid_rowconfigure(8, weight=1)
mainframe_1.grid_rowconfigure(9, weight=1)

ttk.Label(mainframe_1, text="2.Set for PIV.").grid(column=0, row=0, padx=5, pady=5, sticky = W)

def set_aligin():
    Set.need_algin = need_aligin.get() == '1'

need_aligin = StringVar()
ttk.Checkbutton(mainframe_1, text='Need Align', command=set_aligin,variable=need_aligin).grid(column=0, row=1, padx=5, pady=5, sticky = W)
need_aligin.set('1')


def set_norm():
    Set.need_norm = need_norm.get() == '1'

need_norm = StringVar()
ttk.Checkbutton(mainframe_1, text='Need Normalize', command=set_norm, variable=need_norm).grid(column=1, row=1, padx=5, pady=5, sticky = W)
need_norm.set('1')

ttk.Label(mainframe_1, text="Correlation Method").grid(column=0, row=2, padx=5, pady=5, sticky = W)

def set_correlation_method(*args):
    Set.correlation_method = correlation_method.get()

correlation_method = StringVar()
ttk.Combobox(mainframe_1, textvariable=correlation_method, values=('circular', 'linear'), state='readonly').grid(column=1, row=2, padx=5, pady=5, sticky = W)
correlation_method.trace_add("write", set_correlation_method)
correlation_method.set('circular')

ttk.Label(mainframe_1, text="Window Sizes").grid(column=0, row=3, padx=5, pady=5, sticky = W)

def set_windows_size(*args):
    try:
        Set.windowsizes = tuple(np.array(window_size.get().split(), dtype='int'))
    except:
        Set.windowsizes = None

window_size = StringVar()
ttk.Entry(mainframe_1, textvariable=window_size).grid(column=1, row=3, padx=5, pady=5, sticky = W)
window_size.trace_add("write", set_windows_size)
window_size.set('128 64 32')


ttk.Label(mainframe_1, text="Overlaps").grid(column=0, row=4, padx=5, pady=5, sticky = W)

def set_overlap(*args):
    try:
        Set.overlap = tuple(np.array(overlap.get().split(), dtype='int'))
    except:
        Set.overlap = None

overlap = StringVar()
ttk.Entry(mainframe_1, textvariable=overlap).grid(column=1, row=4, padx=5, pady=5, sticky = W)
overlap.trace_add("write", set_overlap)
overlap.set('64 32 16')

ttk.Label(mainframe_1, text="Signal to Noise Threshold").grid(column=0, row=5, padx=5, pady=5, sticky = W)

def set_sig2noise_threshold(*args):
    try:
        Set.sig2noise_threshold = float(sig2noise_threshold.get())
    except:
        Set.sig2noise_threshold = None

sig2noise_threshold = StringVar()
ttk.Entry(mainframe_1, textvariable=sig2noise_threshold).grid(column=1, row=5, padx=5, pady=5, sticky = W)
sig2noise_threshold.trace_add("write", set_sig2noise_threshold)
sig2noise_threshold.set('1.05')

ttk.Label(mainframe_1, text="STD Factor").grid(column=0, row=6, padx=5, pady=5, sticky = W)

def set_std_factor(*args):
    try:
        Set.std_factor = float(std_factor.get())
    except:
        Set.std_factor = None

std_factor = StringVar()
ttk.Entry(mainframe_1, textvariable=std_factor).grid(column=1, row=6, padx=5, pady=5, sticky = W)
std_factor.trace_add("write", set_std_factor)
std_factor.set('10')

ttk.Label(mainframe_1, text="NMT Parameters").grid(column=0, row=7, padx=5, pady=5, sticky = W)

def set_NMT_parameter(*args):
    try:
        a = NMT_parameter.get().split()
        Set.median_threshold = float(a[0])
        Set.median_eps = float(a[1])
        Set.median_size = int(a[2])
    except:
        Set.median_threshold = None

NMT_parameter = StringVar()
ttk.Entry(mainframe_1, textvariable=NMT_parameter).grid(column=1, row=7, padx=5, pady=5, sticky = W)
NMT_parameter.trace_add("write", set_NMT_parameter)
NMT_parameter.set('2.0 0.1 1')

def set_enable_smooth():
    Set.enable_smoothn = enable_smoothn.get() == '1'
enable_smoothn = StringVar()
ttk.Checkbutton(mainframe_1, text='Enable Smooth',command=set_enable_smooth, variable=enable_smoothn).grid(column=0, row=8, padx=5, pady=5, sticky = W)
enable_smoothn.set('1')

def set_smoothn_p(*args):
    try:
        Set.smoothn_p = float(smoothn_p.get())
    except:
        Set.smoothn_p = None
smoothn_p = StringVar()
ttk.Entry(mainframe_1, textvariable=smoothn_p).grid(column=1, row=8, padx=5, pady=5, sticky = W)
smoothn_p.trace_add("write", set_smoothn_p)


ttk.Label(mainframe_1, text="Filter Method").grid(column=0, row=9, padx=5, pady=5, sticky = W)

def set_filter_method(*args):
    Set.filter_method = filter_method.get()

filter_method = StringVar()
ttk.Combobox(mainframe_1, textvariable=filter_method, values=('localmean', 'disk' , 'distance'), state='readonly').grid(column=1, row=9, padx=5, pady=5, sticky = W)
filter_method.trace_add("write", set_filter_method)
filter_method.set('distance')


#----------------------------------------------------------------------------------------------------------#


mainframe_2 = ttk.Frame(root, padding="3 3 12 12")
mainframe_2.grid(column=1, row=0, sticky=(N, W, E, S))
mainframe_2['borderwidth'] = 2
mainframe_2['relief'] = 'raised'

mainframe_2.grid_columnconfigure(0, weight=1)
mainframe_2.grid_columnconfigure(1, weight=1)
mainframe_2.grid_rowconfigure(0, weight=1)
mainframe_2.grid_rowconfigure(1, weight=1)
mainframe_2.grid_rowconfigure(2, weight=1)
mainframe_2.grid_rowconfigure(3, weight=1)
mainframe_2.grid_rowconfigure(4, weight=1)


ttk.Label(mainframe_2, text="3.Set for FTTC.").grid(column=0, row=0, padx=5, pady=5, sticky = W)


ttk.Label(mainframe_2, text="Pixel Size in Micron").grid(column=0, row=1, padx=5, pady=5, sticky = W)

def set_pixelsize(*args):
    try:
        Set.pixelsize = float(pixelsize.get()) * 1e-6
    except:
        Set.pixelsize = None


pixelsize = StringVar()
ttk.Entry(mainframe_2, textvariable=pixelsize).grid(column=1, row=1, padx=5, pady=5, sticky = W)
pixelsize.trace_add("write", set_pixelsize)
pixelsize.set('0.09')

ttk.Label(mainframe_2, text="Young's Modulus in Pascal").grid(column=0, row=2, padx=5, pady=5, sticky = W)

def set_young(*args):
    try:
        Set.young = float(young.get())
    except:
        Set.young = None

young = StringVar()
ttk.Entry(mainframe_2, textvariable=young).grid(column=1, row=2, padx=5, pady=5, sticky = W)
young.trace_add("write", set_young)
young.set('5000')



ttk.Label(mainframe_2, text="Poisson Ratio").grid(column=0, row=3, padx=5, pady=5, sticky = W)

def set_sigma(*args):
    try:
        Set.sigma = float(sigma.get())
    except:
        Set.sigma = None
sigma = StringVar()
ttk.Entry(mainframe_2, textvariable=sigma).grid(column=1, row=3, padx=5, pady=5, sticky = W)
sigma.trace_add("write", set_sigma)
sigma.set('0.5')


ttk.Label(mainframe_2, text="Regularization Parameter").grid(column=0, row=4, padx=5, pady=5, sticky = W)

def set_lamda(*args):
    try:
        Set.lamda = float(lamda.get())
    except:
        Set.lamda = None
lamda = StringVar()
ttk.Entry(mainframe_2, textvariable=lamda).grid(column=1, row=4, padx=5, pady=5, sticky = W)
lamda.trace_add("write", set_lamda)
lamda.set('0.000000000')

#----------------------------------------------------------------------------------------------------------#

mainframe_3 = ttk.Frame(root, padding="3 3 12 12")
mainframe_3.grid(column=1, row=1, sticky=(N, W, E, S))
mainframe_3['borderwidth'] = 2
mainframe_3['relief'] = 'raised'

mainframe_3.grid_columnconfigure(0, weight=1)
mainframe_3.grid_columnconfigure(1, weight=1)
mainframe_3.grid_rowconfigure(0, weight=1)
mainframe_3.grid_rowconfigure(1, weight=1)
mainframe_3.grid_rowconfigure(2, weight=1)
mainframe_3.grid_rowconfigure(3, weight=1)
mainframe_3.grid_rowconfigure(4, weight=1)
mainframe_3.grid_rowconfigure(5, weight=1)
mainframe_3.grid_rowconfigure(6, weight=1)
mainframe_3.grid_rowconfigure(7, weight=1)
mainframe_3.grid_rowconfigure(8, weight=1)
mainframe_3.grid_rowconfigure(9, weight=1)

ttk.Label(mainframe_3, text="4.Set for saving data.").grid(column=0, row=0, padx=5, pady=5, sticky = W)

def set_save_folder():
    folderpath = askdirectory()
    Set.save_folder = folderpath
    save_folder.set(folderpath)
save_folder = StringVar()
ttk.Entry(mainframe_3, textvariable=save_folder, state='readonly' ).grid(row=1, column=0, padx=5, pady=5, sticky = (W,E))
save_folder.set('Folder to save data.')

ttk.Button(mainframe_3, text = 'Open...', command = set_save_folder).grid(row=1, column=1, padx=5, pady=5, sticky = W)

def set_save_PIV():
    Set.save_PIV = save_PIV.get() == '1'

save_PIV = StringVar()
ttk.Checkbutton(mainframe_3, text='Save PIV Plot', command=set_save_PIV,variable=save_PIV).grid(row=2, column=0, padx=5, pady=5, sticky = W)
save_PIV.set('1')


ttk.Label(mainframe_3, text="Plot Height Width").grid(column=0, row=3, padx=5, pady=5, sticky = W)

def set_PIV_height_width(*args):
    try:
        a = PIV_height_width.get().split()
        Set.PIV_height = int(a[0])
        Set.PIV_width = int(a[1])
    except:
        Set.PIV_height = None
        Set.PIV_width = None

PIV_height_width = StringVar()
ttk.Entry(mainframe_3, textvariable=PIV_height_width).grid(column=1, row=3, padx=5, pady=5, sticky = W)
PIV_height_width.trace_add("write", set_PIV_height_width)

ttk.Label(mainframe_3, text="Range of Scale Bar").grid(column=0, row=4, padx=5, pady=5, sticky = W)

def set_PIV_scale_bar(*args):
    try:
        a = PIV_scale_bar.get().split()
        Set.PIV_min_scale = float(a[0])
        Set.PIV_max_scale = float(a[1])
    except:
        Set.PIV_min_scale = -1
        Set.PIV_max_scale = -1

PIV_scale_bar = StringVar()
ttk.Entry(mainframe_3, textvariable=PIV_scale_bar).grid(column=1, row=4, padx=5, pady=5, sticky = W)
PIV_scale_bar.trace_add("write", set_PIV_scale_bar)



ttk.Label(mainframe_3, text="Scale for Length Width Size").grid(column=0, row=5, padx=5, pady=5, sticky = W)

def set_PIV_scale(*args):
    try:
        a = PIV_scale.get().split()
        Set.PIV_scale_for_lenth = float(a[0])
        Set.PIV_scale_of_width = float(a[1])
        Set.PIV_scale_for_arrow = float(a[2])
    except:
        Set.PIV_scale_for_lenth = 1
        Set.PIV_scale_of_width = 1
        Set.PIV_scale_for_arrow = 1

PIV_scale = StringVar()
ttk.Entry(mainframe_3, textvariable=PIV_scale).grid(column=1, row=5, padx=5, pady=5, sticky = W)
PIV_scale.trace_add("write", set_PIV_scale)
PIV_scale.set('1 1 1')


def set_save_FTTC():
    Set.save_FTTC = save_FTTC.get() == '1'

save_FTTC = StringVar()
ttk.Checkbutton(mainframe_3, text='Save FTTC Plot', command=set_save_FTTC,variable=save_FTTC).grid(row=6, column=0, padx=5, pady=5, sticky = W)
save_FTTC.set('1')


ttk.Label(mainframe_3, text="Plot Height Width").grid(column=0, row=7, padx=5, pady=5, sticky = W)

def set_FTTC_height_width(*args):
    try:
        a = FTTC_height_width.get().split()
        Set.FTTC_height = int(a[0])
        Set.FTTC_width = int(a[1])
    except:
        Set.FTTC_height = None
        Set.FTTC_width = None

FTTC_height_width = StringVar()
ttk.Entry(mainframe_3, textvariable=FTTC_height_width).grid(column=1, row=7, padx=5, pady=5, sticky = W)
FTTC_height_width.trace_add("write", set_FTTC_height_width)

ttk.Label(mainframe_3, text="Range of Scale Bar").grid(column=0, row=8, padx=5, pady=5, sticky = W)

def set_FTTC_scale_bar(*args):
    try:
        a = FTTC_scale_bar.get().split()
        Set.FTTC_min_scale = float(a[0])
        Set.FTTC_max_scale = float(a[1])
    except:
        Set.FTTC_min_scale = -1
        Set.FTTC_max_scale = -1
FTTC_scale_bar = StringVar()
ttk.Entry(mainframe_3, textvariable=FTTC_scale_bar).grid(column=1, row=8, padx=5, pady=5, sticky = W)
FTTC_scale_bar.trace_add("write", set_FTTC_scale_bar)


ttk.Label(mainframe_3, text="Scale for Length Width Size").grid(column=0, row=9, padx=5, pady=5, sticky = W)

def set_FTTC_scale(*args):
    try:
        a = FTTC_scale.get().split()
        Set.FTTC_scale_for_lenth = float(a[0])
        Set.FTTC_scale_of_width = float(a[1])
        Set.FTTC_scale_for_arrow = float(a[2])
    except:
        Set.FTTC_scale_for_lenth = 1
        Set.FTTC_scale_of_width = 1
        Set.FTTC_scale_for_arrow = 1
FTTC_scale = StringVar()
ttk.Entry(mainframe_3, textvariable=FTTC_scale).grid(column=1, row=9, padx=5, pady=5, sticky = W)
FTTC_scale.trace_add("write", set_FTTC_scale)
FTTC_scale.set('1 1 1')

#----------------------------------------------------------------------------------------------------------#
mainframe_4 = ttk.Frame(root, padding="3 3 12 12")
mainframe_4.grid(column=0, row=3, columnspan=2, sticky=(N, W, E, S))
mainframe_4['borderwidth'] = 2
mainframe_4['relief'] = 'raised'
mainframe_4.grid_columnconfigure(0, weight=1)

pb = StringVar()
ttk.Progressbar(mainframe_4, value=0, mode="determinate", variable=pb).grid(column=0, row=0, padx=5, pady=5, sticky=(W,E))
pb.set(0)

strat_button = ttk.Button(mainframe_4, text='Start', command = check_parameters)
strat_button.grid(column=1, row=0, padx=5, pady=5, sticky = E)
#----------------------------------------------------------------------------------------------------------#
root.resizable(0,0)


if __name__ == '__main__':
    root.mainloop()