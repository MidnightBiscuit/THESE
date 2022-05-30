import numpy as np
from matplotlib import pylab, mlab, pyplot
from pylab import *
from IPython.core.pylabtools import figsize, getfigs
plt = pyplot

import os
from os import listdir
from os.path import isfile, join
import sys

from pathlib import Path
import re
from natsort import natsorted
import timeit

# sys.path.append('/home/adrian/Documents/Programmes/Fortran_Jofre')
# from function_jofre import energy_lost

# sélection des fichiers (Python 3)
import tkinter as tk
from tkinter import filedialog

def load_file_GUI(dir_string):

    root = tk.Tk()
    root.withdraw()

    file_path = filedialog.askopenfilename(initialdir = dir_string,
                                       multiple=True)
    return file_path

def data_retrieve(all_subdir,points_and_coord, condition_parameters, slash_cfg,**kwargs):

    myslashpos = slash_cfg[0]
    slashcond = slash_cfg[1]

    # determining number of elements on each repetition
    try:
        num_runs = ['Try'+str(kwargs.get('forcenumtry'))]
        num_runs_aux = [runs[myslashpos[slashcond+1]+1:] for runs in all_subdir if list(points_and_coord.keys())[0] in runs]
        num_runs_aux = list(dict.fromkeys(num_runs_aux))
        num_runs_aux = len(num_runs_aux)
        all_subdir = [x for w,x in enumerate(all_subdir) if (w-kwargs.get('forcenumtry'))%num_runs_aux==0]
    except:
        num_runs = [runs[myslashpos[slashcond+1]+1:] for runs in all_subdir if list(points_and_coord.keys())[0] in runs]
        num_runs = list(dict.fromkeys(num_runs))

    # number of repetitions
    print('> Points |',len(points_and_coord))
    print('> Simulations pour chaque point |', num_runs)
    
    data0 = [[] for i in range(len(points_and_coord))] # Size
    data1 = [[] for i in range(len(points_and_coord))] # Temperature
    data2 = [[] for i in range(len(points_and_coord))] # xva
    data3 = [[] for i in range(len(points_and_coord))] # trj
    data_address = [[] for i in range(len(points_and_coord))]

    # Variables à deux coordonnées : [point, try]
    shapevar = (len(points_and_coord),len(num_runs))
    Time        = []
    Temperature = []
    Size_w      = []
    Size2       = []    
    t0 = time.clock()
    print("Hello")

    # write variables
    for k, address in enumerate(all_subdir):

    # in-loop variables
        pnt = k // len(num_runs)  # actual point
        rep = k  % len(num_runs)          # actual repetition

        # get only .dat files in each simulation directory
        print(address)
        onlyfiles = [f for f in listdir(address) if isfile(join(address, f)) and not "xva" in f and ".dat" in f]
        print(onlyfiles)
        # build path file
        data0[pnt].append('{}/{}'.format(address,sort(onlyfiles)[0].strip('.dat')))
        data1[pnt].append('{}/{}'.format(address,sort(onlyfiles)[1].strip('.dat')))
        data_address[pnt].append(address)

        # load fluorescence and T
        Time.append( loadtxt(str(data1[pnt][0])+'.dat',unpack=True)[0] )
        Temperature.append( loadtxt(str(data1[pnt][0])+'.dat',unpack=True)[4:7] )
        Size_w.append( loadtxt(str(data0[pnt][0])+'.dat',unpack=True)[0:3] )
        Size2.append( loadtxt(str(data0[pnt][0])+'.dat',unpack=True)[3] )
        
        onlyfiles = [f for f in listdir(address) if isfile(join(address, f)) and "xva" in f]
        data2[pnt].append('{}/{}'.format(address,sort(onlyfiles)[0].strip('.bin')))
        
        onlyfiles = [f for f in listdir(address) if isfile(join(address, f)) and "trj" in f]
        data3[pnt].append('{}/{}'.format(address,sort(onlyfiles)[0].strip('.bin')))
        
        # load cloud size before injection
        if not(rep % len(num_runs)):
            print( "Point n°", pnt )
        
        print(f'{pnt:02}','-',f'{rep:02}',' > ',data0[pnt][rep])
        
        if k == len(points_and_coord)*len(num_runs) - 1:
            break

    t1 = time.clock() - t0
    print("Time elapsed: ", t1, 's') # CPU seconds elapsed (floating point)
    print("Time elapsed: ", t1/60, 'm') # CPU seconds elapsed (floating point)
    
    data_name = [data_address, data0, data1, data2, data3]
    
#    print('{}/{}'.format(address,sort(onlyfiles)[0].strip('.dat')))
#    print('{}/{}'.format(address,sort(onlyfiles)[1].strip('.dat')))
    
    return data_name, num_runs, Time, Temperature, Size_w, Size2

def load_gui(filter_nocomplete):

    # sélection des fichiers (Python 3)
    # SELECTIONNER UN FICHIER DE DONNÉES AU FOND DE L'ARBORESCENCE
    # Temp_SimuType0_N01024_Vrf0064_Udc0.5000D+00V_D1.0_S1.0RFG.dat

#     root = tk.Tk()
#     root.withdraw()

#     # Sélectionner un fichier du dossier pour récupérer le répertoire et les noms
#     file_path = filedialog.askopenfilename(initialdir = '/home/adrian/Documents/Programmes/Fortran_Jofre',
#                                            multiple=None)
    
    file_path = load_file_GUI('/home/adrian/Documents')[0]
    
    dir_path = Path(file_path).parent.absolute()
    work_rep = file_path[:len(str(dir_path))]
    filename = file_path[len(str(dir_path))+1:]

    print('{}{}'.format('> Répertoire : ',work_rep))
    print('{}{}'.format('> Filename : ',filename))
    
    ## Recover all subdirectories hosting simulation data

    # find all '/' in path
    myslashpos = [m.start() for m in re.finditer('/', work_rep)]
    # print(myslashpos[6]) # choose the right level

    if 'home' in file_path:
        if 'Hobitton' or 'Rivendel' in file_path :
            slashcond = -1 # the slash number just before runs (after date)
        else:
            slashcond = -1 # the slash number just before runs (after date)
    else:
        slashcond = -1

    print('> myslashpos |',myslashpos)
    print('> slashcond |',slashcond)

    # Getting a naturally sorted list of all subdirectories in some directory
    # throw all the intermediary subdirectories that are not hosting data
    all_subdir = [x[0] for x in os.walk(work_rep[0:myslashpos[slashcond]])]
    all_subdir = natsorted(all_subdir)

    # keep only the long path (by remove short ones)
    all_subdir = [i for i in all_subdir if len(i) >= len(str(dir_path))-1]

    if filter_nocomplete == 1:
        # Name of all points (conditions) 
        all_points = [point[myslashpos[slashcond] + 1:myslashpos[slashcond + 1]] for point in all_subdir]
        all_points = list(dict.fromkeys(all_points))
        # how much repetition for first point : usually the first point is complete
        # i'm not going to plot the data while not even the first point is complete ...
        num_rep = len([k for k in all_subdir if all_points[0] in k])
        
        # determine which points are not complete
        # the way the code is done usually leads to
        # only one uncomplete point but in any case
        # I handle the possibility for multiple uncomplete.
        all_subdir_temp = []
        pts_to_delete = []
        deleted_points = 0
        
        # sweep all points
        for i,j in enumerate(all_points):
            temp = [k for k in all_subdir if j in k] # all addresses given one point condition
            if len(temp) != num_rep: # compare number of repetition for each point
                                     # if different from the reference (num_rep)
                                     # I will remove this point
                pts_to_delete.append(j)
                
        # removing uncomplete points from the list all_subdir
        if len(pts_to_delete) != 0:
            for i,w in enumerate(all_subdir): # sweep all points
                # ~ print(w)
                for j,x in enumerate(pts_to_delete): # sweep all delete conditions
                    # ~ print(x)
                    if x not in w: # if the current point is not complete, delete from all_subdir
                        all_subdir_temp.append(w)
                    else:
                        deleted_points += 1
                        # print('===== ^ REMOVED ^ =====')
            all_subdir = all_subdir_temp
        del all_subdir_temp
                    
        print('Points deleted because they were not complete',pts_to_delete,'  '+str(deleted_points)+' pt(s)')
        print('Total number of data directories',len(all_subdir))
    else:
        print('No points deleted because they were not complete')
        print('Total number of data directories',len(all_subdir))        

    # all_subdir ce sont tous les répertoires contenant des donnés à analyser
    
    file_cfg = [file_path, dir_path, work_rep, filename]
    slash_cfg = [myslashpos, slashcond]
        
    return file_cfg, slash_cfg, all_subdir
    
def simu_conditions(all_subdir, myslashpos, slashcond,filename):

    # All points of simulation
    all_points = [point[myslashpos[slashcond] + 1:myslashpos[slashcond + 1]] for point in all_subdir]
    all_points = list(dict.fromkeys(all_points))
    print(all_points)

    # Name of the conditions
    condition_separator_position = [m.start() for m in re.finditer('_', all_points[0])]
    if len(condition_separator_position) == 0:
        condition_name = [re.sub('[0-9]+','',all_points[0])]
        condition_numb = [re.sub('[a-zA-Z]','',kk) for _,kk in enumerate(all_points)]
        points_and_coord = [k for j,k in enumerate(all_points)]
        points_and_coord = {y: condition_numb[x] for x,y in enumerate(all_points)}

    else:
        condition_name = [[] for k in range(len(condition_separator_position) + 1)]
        for k, m in enumerate(condition_separator_position):
            condition_name[k] = re.sub('[0-9]+', '', all_points[0][m - condition_separator_position[0]:m])
        condition_name[-1] = re.sub('[0-9]+', '', all_points[0][m + 1:])
        
        # Put together points with their coordinates        
        points_and_coord = {}
        for k,pt in enumerate(all_points):
            print(f'{k:03.0f}','>',pt)
            w = condition_separator_position[0]
            temp = []
            for i,j in enumerate(condition_name):
                temp_cond_pos = pt.find(j)
                try:
                    temp_sep_pos = condition_separator_position[i]
                    temp_cond_num = re.sub(j,'',pt[temp_cond_pos:temp_sep_pos])
                except:
                    temp_cond_num = re.sub(j,'',pt[temp_cond_pos:])
                temp.append(temp_cond_num)
            points_and_coord[pt] = temp
        # {'DC01_RF08': ['01', '08'], 'DC01_RF09': ['01', '09'], 'DC01_RF10': ['01', '10'], ... }

    print('> condition names',condition_name)
    print('> number of points',len(all_points))

    # Conditions de la simulation
    nions_pos = filename.find('_N')            # position nombre d'ions dans filename
    N = int(filename[nions_pos+2:nions_pos+7]) # nombre d'ions
    e_GMol = 50   # energie GMol en eV

    print(f'> N_ions = {N}')
    print(f'> e_GMol = {e_GMol}')

    condition_parameters = [condition_name, N, e_GMol]

    return points_and_coord, condition_parameters

def load_xva_init_bin_Ver(str_load, N_ions):      

    fid = open(str_load+'.bin',  'rb')   
    
    e    = fromfile(fid,  int32,   1   )        # Read record start tag   
    iRF  = fromfile(fid,  int32, 3)   
    junk = fromfile(fid,  int32,   1   )        # Read record stop tag   

    f    = fromfile(fid,  int32,   1   )        # Read record start tag 
    save_xva = fromfile(fid,  float64, 12*N_ions)   
    junk = fromfile(fid,  int32,   1   )        # Read record stop tag

    fid.close   
    
    # print('len save_T', junk)   

    save_xva = reshape(save_xva, (12,N_ions), order='F')   
#     print(save_xva[0])   

    return iRF, save_xva

def beta_continue_alamano(a,q,beta_guess): # a l'air ok
    C1 = q**2/( (beta_guess+2)**2 - a - q**2/( (beta_guess+4)**2 - a - q**2/( (beta_guess+6)**2 - a - q**2/( (beta_guess+8)**2 - a - q**2/( (beta_guess+10)**2 - a ) ) ) ) )
    C2 = q**2/( (beta_guess-2)**2 - a - q**2/( (beta_guess-4)**2 - a - q**2/( (beta_guess-6)**2 - a - q**2/( (beta_guess-8)**2 - a - q**2/( (beta_guess-10)**2 - a  ) ) ) ) )
    return a + C1 + C2
# ======================================================================== #
# pour les couleurs
from matplotlib import colors as mcolors
def cc(arg,alpha):
    '''
    Shorthand to convert 'named' colors to rgba format at x% opacity.
    '''
    return mcolors.to_rgba(arg, alpha=alpha)

# convert H:M:s time to sec
def get_sec(time_str):
    h, m, s = time_str.split(':')
#     return h,m,s
    return int(h) * 3600 + int(m) * 60 + int(s)
