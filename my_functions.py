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
# from natsort import natsorted
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

def size_data(file_path,row_skip,col_to_read,delim):
    all_size = []
    m = 0
    for k in range(0, len(file_path)):
        temp_size = loadtxt(file_path[k], delimiter=delim,
                            skiprows=row_skip, usecols=col_to_read,
                            unpack=True)[0].size
        all_size.append(temp_size)
        if temp_size > m:
            m = temp_size
    return all_size,m

# import données    /usr/share/myspell/dicts/en_GB.dic
def import_data(file_path,row_skip,col_to_read,delim):

    dico = {}                       # creation dictionnaire
    for i in file_path:              # pour associer à chaque fichier
        dico[i] = []                # un nombre de colonnes dépendant de col_to_read
# dico de la forme
# dico = {'my_file1.csv': array([1,2,3]),
#          'my_file2.csv': array([2,4,6]),
#          'my_file3.csv': array([5,10,15])}
# récup valeur dans dico : dico.get('my_file1.csv')

    for k in range(0, len(file_path)):
        dico[file_path[k]] = loadtxt(file_path[k], delimiter=delim,encoding ='utf-8',
                                      skiprows=row_skip, usecols=col_to_read,
                                      unpack=True)
    return dico

def convert_dico_to_var(dico):

    file_path = list(dico.keys()) # all the keys from dico

    # création variables courantes
    n = len(dico)  # nombre de fichiers
    m = len(dico.get(file_path[0])[0])  # longueur colones
    p = len(dico.get(file_path[0])[1:])  # nombre de col : temps + canaux oscillo
    shape = (n, m)
    TP = zeros(shape)  # time
    CH = zeros((n, p, m))  # les canaux de l'oscillo (p-1)  CH[file:channel:time]

    for k in range(0, n):
        TP[k, :] = dico.get(file_path[k])[0, :]
        for l in range(1, p+1):
            CH[k, l-1, :] = dico.get(file_path[k])[l, :]

    return TP,CH

##### FUNCTION_JOFRE.IPY #####

def load_T_and_PM_simu(str_load):
    data = loadtxt('{}.dat'.format(str_load),comments='%');
    return data[:,0],data[:,1:4], data[:,4:7], data[:,7] # tt, T_CM, T_aux, PM
    
def load_T_and_PM_simu_noGMOL(str_load):
    data = loadtxt('{}.dat'.format(str_load),comments='%');
    return data[:,0],data[:,1:4] # tt, T_CM, T_aux, PM

def load_xyz_init_bin_DP(str_load):
    aux_info = loadtxt(str_load+'.info',comments='%');
    n_ions1 = int(aux_info[0])
    n_ions2 = int(aux_info[1])
    n_ions = n_ions1 + n_ions2
#    print(n_ions1, n_ions2)
    imax = 12 * n_ions

#    print('extra_time', aux_info[19])
    
    fid = open(str_load+'.bin', 'rb')
    junk = fromfile(fid, int32,1)        # Read record start tag
    aux  = fromfile(fid, int32, 3);
    junk = fromfile(fid, int32,1)        # Read record stop tag

    junk = fromfile(fid, int32,1)        # Read record start tag
    aux  = fromfile(fid, float64, imax);
    junk = fromfile(fid, int32,1)        # Read record stop tag

    fid.close
    aux = reshape(aux, (12, n_ions),order='F')
    r_LC = 1.e3*aux[0:3,0:n_ions1]
    v_LC = aux[3:6,0:n_ions1]

    #xyz = 1.e3*aux[0:3,:]
    a_LC = aux[6:9,0:n_ions1]
    v_LC_avg = aux[9:12,0:n_ions1]

    return r_LC,v_LC,a_LC,v_LC_avg

def plot_XYZ(file_name,fig_name='2',fig_title='XYZ'):
    r_LC,v_LC,a_LC,v_LC_avg = load_xyz_init_bin_DP(file_name)
    figure(fig_name); clf()
    title(fig_title)
    subplot(211,aspect=1.0)
    plot(r_LC [0,:]*1e3,r_LC [1,:]*1e3,'8',color='xkcd:purplish blue')
    xlabel('x[µm]')
    ylabel('y[µm]')
    grid()

    # subplot(212,aspect=1.0)
    subplot(212)
    plot(r_LC [2,:]*1e3,r_LC [0,:]*1e3,'8',color='xkcd:purplish blue')
    xlabel('z[µm]')
    ylabel('x[µm]')
    grid()
    
    tight_layout()
    
# def plot_T_and_PM_Init_Inje_Evol(file_dir,file1, file2, file3,fig_name='3'):

def plot_T_and_PM_Init_Inje_Evol(file_dir2,file_name,flag_plot,fig_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.5*1e-3,2e4))
    ylim2 = kwargs.get('ylim2', (-2,50))
    
    i_aux = file_name.find('_N')
    file0 = file_dir2+'Temp_3D_Harmo_N1024_T500uK_F0.15D-18Kg_s_5'
    file1 = 'SimuTypeQ'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]
    file2 = 'SimuType4_01' + file_name[i_aux:]
    file3 = 'SimuType2_01' + file_name[i_aux:]

    N_ions, j_save, dt_j_save_next, eta, Temp, save_T = load_Temp_init_bin_Lan(file0, flag_print=0)
    tt0 = save_T[:,0]
    T_CM0 = save_T[:,1:4]
    T_aux0 = save_T[:,4:7]
    # tt0, T_CM0, T_aux0, PM0 = load_T_and_PM_simu(file_dir2+file0)
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
    tt2, T_CM2, T_aux2, PM2 = load_T_and_PM_simu(file_dir2+'Temp_'+file2+'50eV')
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'Temp_'+file3+'50eV')

    aux = mean(PM1[-100:])
    PM_variation = ( aux - mean(PM3[-100:]) ) / aux
    T_variation  = mean(T_aux3[-100:,0]) + mean(T_aux3[-100:,1]) + mean(T_aux3[-100:,2])
            

    # Auxiliary arrays:
    t_aux1 = array([tt2[ 0],tt2[ 0]])
    t_aux2 = array([tt2[-1],tt2[-1]])
    y1_aux = array([1.0e-3 ,1.0e-1 ])
#     y2_aux = array([0 ,20 ])
    y2_aux = array([0 ,50 ])

    tt    = concatenate( ( tt0,   tt1,   tt2,   tt3) )
    tt_laser = concatenate( ( tt1,   tt2,   tt3) )
    T_CM  = concatenate( (T_CM0,T_CM1,T_CM2,T_CM3) )
    T_aux = concatenate( (T_aux0,T_aux1,T_aux2,T_aux3) )
    PM    = concatenate( (PM1,PM2,PM3) )
    SNR = np.abs( aux - mean(PM3[-100:]) )/np.sqrt(aux) # Sig-to-Noi ratio considering Poisson noise
    
    if flag_plot == 1 :
        #fig_name = file_name[-9:]
        fig = figure(fig_name); clf()
        ax1 = subplot(211)
        semilogy(tt*1.e3,T_aux[:,0], label='Tx')
        semilogy(tt*1.e3,T_aux[:,1], label='Ty')
        semilogy(tt*1.e3,T_aux[:,2], label='Tz')
        semilogy(t_aux1*1.e3,y1_aux,'r')
        semilogy(t_aux2*1.e3,y1_aux,'r')
        ax1.grid()

        legend()
        # ~ xlabel('time[ms]')
        # ~ ylabel('T[K]')
        plt.setp(ax1.get_xticklabels(),visible=False)

        ax2 = subplot(212,sharex=ax1)
        plot(tt_laser*1.e3,PM[:])
        plot(t_aux1*1.e3,y2_aux,'r')
        plot(t_aux2*1.e3,y2_aux,'r')
        ax2.grid()
        
        xlabel('time[ms]')
        ylabel('Counts')
                               
        ax1.set_xlim(xlim1)
        ax1.set_ylim(ylim1)
        ax2.set_ylim(ylim2)
        plt.tight_layout()
        subplots_adjust(hspace=0.015)
        
    temps = [tt,t_aux1, t_aux2]
    fluo = [PM, PM_variation, SNR]
    temperature = [T_aux, T_variation, T_CM]
        
    return temps, temperature, fluo
    
def plot_T_and_PM_InitQ_Inje_Evol(file_dir2,file_name,flag_plot,fig_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.5*1e-3,2e4))
    ylim2 = kwargs.get('ylim2', (-2,50))
    
    i_aux = file_name.find('_N')
    # file0 = 'Temp_3D_Harmo'+ file_name[i_aux:]
    file1 = 'SimuTypeQ'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]
    file2 = 'SimuType4_01' + file_name[i_aux:]
    file3 = 'SimuType2_01' + file_name[i_aux:]

    # tt0, T_CM0, T_aux0, PM0 = load_T_and_PM_simu(file_dir2+file0)
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
    tt2, T_CM2, T_aux2, PM2 = load_T_and_PM_simu(file_dir2+'Temp_'+file2+'50eV')
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'Temp_'+file3+'50eV')

    aux = mean(PM1[-100:])
    PM_variation = ( aux - mean(PM3[-100:]) ) / aux
    T_variation  = mean(T_aux3[-100:,0]) + mean(T_aux3[-100:,1]) + mean(T_aux3[-100:,2])
            

    # Auxiliary arrays:
    t_aux1 = array([tt2[ 0],tt2[ 0]])
    t_aux2 = array([tt2[-1],tt2[-1]])
    y1_aux = array([1.0e-3 ,1.0e-1 ])
#     y2_aux = array([0 ,20 ])
    y2_aux = array([0 ,50 ])

    tt    = concatenate( (   tt1,   tt2,   tt3) )
    T_CM  = concatenate( (T_CM1,T_CM2,T_CM3) )
    T_aux = concatenate( (T_aux1,T_aux2,T_aux3) )
    PM    = concatenate( (PM1,PM2,PM3) )
    SNR = np.abs( aux - mean(PM3[-100:]) )/np.sqrt(aux) # Sig-to-Noi ratio considering Poisson noise
    
    if flag_plot == 1 :
        #fig_name = file_name[-9:]
        fig = figure(fig_name); clf()
        ax1 = subplot(211)
        semilogy(tt*1.e3,T_aux[:,0], label='Tx')
        semilogy(tt*1.e3,T_aux[:,1], label='Ty')
        semilogy(tt*1.e3,T_aux[:,2], label='Tz')
        semilogy(t_aux1*1.e3,y1_aux,'r')
        semilogy(t_aux2*1.e3,y1_aux,'r')
        ax1.grid()

        legend()
        # ~ xlabel('time[ms]')
        # ~ ylabel('T[K]')
        plt.setp(ax1.get_xticklabels(),visible=False)

        ax2 = subplot(212,sharex=ax1)
        plot(tt*1.e3,PM[:])
        plot(t_aux1*1.e3,y2_aux,'r')
        plot(t_aux2*1.e3,y2_aux,'r')
        ax2.grid()
        
        xlabel('time[ms]')
        ylabel('Counts')
                               
        ax1.set_xlim(xlim1)
        ax1.set_ylim(ylim1)
        ax2.set_ylim(ylim2)
        plt.tight_layout()
        subplots_adjust(hspace=0.015)
        
    temps = [tt,t_aux1, t_aux2]
    fluo = [PM, PM_variation, SNR]
    temperature = [T_aux, T_variation, T_CM]
        
    return temps, temperature, fluo

def plot_T_and_PM_InitQ_Inje(file_dir2,file_name,flag_plot,fig_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.5*1e-3,2e4))
    ylim2 = kwargs.get('ylim2', (-2,50))
    
    i_aux = file_name.find('_N')
    # file0 = 'Temp_3D_Harmo'+ file_name[i_aux:]
    file1 = 'SimuTypeQ'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]
    file2 = 'SimuType4_01' + file_name[i_aux:]

    # tt0, T_CM0, T_aux0, PM0 = load_T_and_PM_simu(file_dir2+file0)
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
    tt2, T_CM2, T_aux2, PM2 = load_T_and_PM_simu(file_dir2+'Temp_'+file2+'50eV')

    aux = mean(PM1[-100:])            

    # Auxiliary arrays:
    t_aux1 = array([tt2[ 0],tt2[ 0]])
    t_aux2 = array([tt2[-1],tt2[-1]])
    y1_aux = array([1.0e-3 ,1.0e-1 ])
#     y2_aux = array([0 ,20 ])
    y2_aux = array([0 ,50 ])

    tt    = concatenate( (   tt1,   tt2) )
    T_CM  = concatenate( (T_CM1,T_CM2) )
    T_aux = concatenate( (T_aux1,T_aux2) )
    PM    = concatenate( (PM1,PM2) )
    
    if flag_plot == 1 :
        #fig_name = file_name[-9:]
        fig = figure(fig_name); clf()
        ax1 = subplot(211)
        semilogy(tt*1.e3,T_aux[:,0], label='Tx')
        semilogy(tt*1.e3,T_aux[:,1], label='Ty')
        semilogy(tt*1.e3,T_aux[:,2], label='Tz')
        semilogy(t_aux1*1.e3,y1_aux,'r')
        semilogy(t_aux2*1.e3,y1_aux,'r')
        ax1.grid()

        legend()
        # ~ xlabel('time[ms]')
        # ~ ylabel('T[K]')
        plt.setp(ax1.get_xticklabels(),visible=False)

        ax2 = subplot(212,sharex=ax1)
        plot(tt*1.e3,PM[:])
        plot(t_aux1*1.e3,y2_aux,'r')
        plot(t_aux2*1.e3,y2_aux,'r')
        ax2.grid()
        
        xlabel('time[ms]')
        ylabel('Counts')
                               
        ax1.set_xlim(xlim1)
        ax1.set_ylim(ylim1)
        ax2.set_ylim(ylim2)
        plt.tight_layout()
        subplots_adjust(hspace=0.015)
        
    temps = [tt,t_aux1, t_aux2]
    fluo = [PM, 0, 0]
    temperature = [T_aux, T_variation, T_CM]
        
    return temps, temperature, fluo


def plot_T_and_PM_Init(file_dir2,file_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.5*1e-3,2e4))
    ylim2 = kwargs.get('ylim2', (-2,50))
    
    i_aux = file_name.find('_N')
    # file0 = 'Temp_3D_Harmo'+ file_name[i_aux:]
    file1 = 'SimuType0'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]

    # tt0, T_CM0, T_aux0, PM0 = load_T_and_PM_simu(file_dir2+file0)
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
            

    # Auxiliary arrays:
    t_aux1 = array([tt1[-1],tt1[-1]])
    y1_aux = array([1.0e-3 ,1.0e-1 ])
        
    return tt1, T_CM1, T_aux1, PM1

def plot_T_and_PM_Quen(file_dir2,file_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.5*1e-3,2e4))
    ylim2 = kwargs.get('ylim2', (-2,50))
    
    i_aux = file_name.find('_N')
    # file0 = 'Temp_3D_Harmo'+ file_name[i_aux:]
    file1 = 'SimuTypeQ'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]

    # tt0, T_CM0, T_aux0, PM0 = load_T_and_PM_simu(file_dir2+file0)
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
            

    # Auxiliary arrays:
    t_aux1 = array([tt1[-1],tt1[-1]])
    y1_aux = array([1.0e-3 ,1.0e-1 ])
        
    return tt1, T_CM1, T_aux1, PM1

def plot_T_and_PM_Inje(file_dir2,file_name,**kwargs):
   
    i_aux = file_name.find('_N')
    file2 = 'SimuType4_01' + file_name[i_aux:]

    tt2, T_CM2, T_aux2, PM2 = load_T_and_PM_simu(file_dir2+'/Temp_'+file2+'50eV')

    T_variation  = mean(T_aux2[-25:,0]) + mean(T_aux2[-25:,1]) + mean(T_aux2[-25:,2])

    # Auxiliary arrays:
    t_aux1 = array([tt2[ 0],tt2[ 0]])
    t_aux2 = array([tt2[-1],tt2[-1]])
    y1_aux = array([1.0e-3 ,1.0e-1 ])
#     y2_aux = array([0 ,20 ])
    y2_aux = array([0 ,50 ])
    
    return tt2, T_CM2, T_aux2, PM2

def plot_T_and_PM_Evol(file_dir2,file_name,**kwargs):

    i_aux = file_name.find('_N')
    file3 = 'SimuType2_01' + file_name[i_aux:]

    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'Temp_'+file3+'50eV')

    T_variation  = mean(T_aux3[-100:,0]) + mean(T_aux3[-100:,1]) + mean(T_aux3[-100:,2])
        
    return tt3, T_CM3, T_aux3, PM3

def plot_T_and_PM_Init_RFrelax_AfterCooling_Evol(file_dir2,file_name,flag_plot,fig_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.5*1e-3,2e4))
    ylim2 = kwargs.get('ylim2', (-2,50))
    
    i_aux = file_name.find('_N')
    file1 = 'SimuType0'    + file_name[i_aux:]
    file2 = 'SimuType4_01' + file_name[i_aux:]
    file3 = 'SimuType2_01' + file_name[i_aux:]
    file4 = 'SimuType6_01' + file_name[i_aux:]

    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
    tt2, T_CM2, T_aux2, PM2 = load_T_and_PM_simu(file_dir2+'Temp_'+file2+'50eV')
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'Temp_'+file3+'50eV')
    try:
        tt4, T_CM4, T_aux4, PM4 = load_T_and_PM_simu(file_dir2+'Temp_'+file4+'50eV')        
    except:
        tt4, T_CM4, T_aux4, PM4 = load_T_and_PM_simu(file_dir2+'Temp_'+file4)

    T_variation  = mean(T_aux2[-100:,0])+mean(T_aux2[-100:,1])+mean(T_aux2[-100:,2]) - mean(T_aux1[-100:,0])+mean(T_aux1[-100:,1])+mean(T_aux1[-100:,2])
    aux = mean(PM1[-100:])
    PM_variation = ( aux - mean(PM3[-100:]) ) / aux
    SNR = np.abs( aux - mean(PM3[-100:]) )/np.sqrt(aux) # Sig-to-Noi ratio considering Poisson noise

    # Auxiliary arrays:
    t_aux1 = array([tt4[ 0],tt4[ 0]])
    t_aux2 = array([tt4[-1],tt4[-1]])
    y1_aux = array([1.0e-4 ,3.0e4 ])
#     y2_aux = array([0 ,20 ])
    y2_aux = array([0 ,50 ])

    tt    = concatenate( (   tt1,   tt4) )
    T_CM  = concatenate( (T_CM1, T_CM4) )
    T_aux = concatenate( (T_aux1,T_aux4) )
    PM    = concatenate( (PM1,PM4) )
    
    if flag_plot == 1 :
        #fig_name = file_name[-9:]
        fig = figure(fig_name); clf()
        ax1 = subplot(111)
        semilogy(tt*1.e3,T_aux[:,0], label='Tx')
        semilogy(tt*1.e3,T_aux[:,1], label='Ty')
        semilogy(tt*1.e3,T_aux[:,2], label='Tz')
        semilogy(t_aux1*1.e3,y1_aux,'r')
        semilogy(t_aux2*1.e3,y1_aux,'r')
        ax1.grid()
        # ax1.set_xticks([1e-4,1e-2,1,1e2,1e4])

        legend()
        xlabel('time[ms]')
        ylabel('T[K]')
        
        fig.set_size_inches(11.69, 8.27)
        plt.tight_layout()
        
    temps = [tt,t_aux1, t_aux2]
    fluo = [PM, PM_variation, SNR]
    temperature = [T_aux, T_variation, T_CM]
        
    return temps, temperature, fluo

def plot_T_and_PM_Init_RFrelax_AfterInj_Evol(file_dir2,file_name,flag_plot,fig_name,**kwargs):
    
    # ~ xlim1 = (-0.1,6)
    # ~ ylim1 = (0.5*1e-3,5e3)
    # ~ ylim2 = (-2,120)
            
    xlim1 = kwargs.get('xlim1', (-0.1,6))
    ylim1 = kwargs.get('ylim1', (0.2*1e-3,2e4))
    ylim3 = kwargs.get('ylim2', (-2,85))
    
    i_aux = file_name.find('_N')
    file1 = 'SimuType0'    + file_name[i_aux:] ########## file1 = 'SimuType0'    + file_name[i_aux:17+36]
    file2 = 'SimuType4_01' + file_name[i_aux:]
    file3 = 'SimuType2_01' + file_name[i_aux:]
    file4 = 'SimuType6_01' + file_name[i_aux:]

    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'Temp_'+file1)
    tt2, T_CM2, T_aux2, PM2 = load_T_and_PM_simu(file_dir2+'Temp_'+file2+'50eV')
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'Temp_'+file3+'50eV')
    try:
        tt4, T_CM4, T_aux4, PM4 = load_T_and_PM_simu(file_dir2+'Temp_'+file4+'50eV')        
    except:
        tt4, T_CM4, T_aux4, PM4 = load_T_and_PM_simu(file_dir2+'Temp_'+file4)
    
    aux = mean(PM1[-100:])
    PM_variation = ( aux - mean(PM3[-100:]) ) / aux
    T_variation  = mean(T_aux3[-100:,0]) + mean(T_aux3[-100:,1]) + mean(T_aux3[-100:,2])
    SNR = np.abs( aux - mean(PM3[-100:]) )/np.sqrt(aux) # Sig-to-Noi ratio considering Poisson noise

    # Auxiliary arrays:
    t_aux1 = array([tt2[ 0],tt2[ 0]])
    t_aux2 = array([tt2[-1],tt2[-1]])
    y1_aux = array([1.0e-3 ,1.0 ])
#     y2_aux = array([0 ,20 ])
    y2_aux = array([0 ,50 ])

    tt    = concatenate( (   tt1,   tt2,   tt3) )
    T_aux = concatenate( (T_aux1,T_aux2,T_aux3) )
    T_CM  = concatenate( (T_CM1,T_CM2,T_CM3) )
    PM    = concatenate( (PM1,PM2,PM3) )
    
    tt_relax    = concatenate( (   tt1,   tt2,   tt4) )
    T_aux_relax = concatenate( (T_aux1,T_aux2,T_aux4) )
    PM_relax    = concatenate( (PM1,PM2,PM4) )
    
    if flag_plot == 1 :
        #fig_name = file_name[-9:]
        fig = figure(fig_name); clf()
        ax1 = subplot(311)
        semilogy(tt*1.e3,T_aux[:,0], label='Tx')
        semilogy(tt*1.e3,T_aux[:,1], label='Ty')
        semilogy(tt*1.e3,T_aux[:,2], label='Tz')
        semilogy(t_aux1*1.e3,y1_aux,'r')
        semilogy(t_aux2*1.e3,y1_aux,'r')
        ax1.set_yticks([1e-4,1e-2,1,1e2,1e4])
        ax1.grid()
        # annotate('Laser ON', xy=(0.5,350), xycoords='data',
            # size=24, ha='left', va='top', color='xkcd:azul',
            # bbox=dict(boxstyle='round', fc='white',edgecolor='xkcd:azul'))
        

        legend(title='Laser ON',fontsize=18)
        # ~ xlabel('time[ms]')
        # ~ ylabel('T[K]')
        plt.setp(ax1.get_xticklabels(),visible=False)

        ax2 = subplot(312,sharex=ax1,sharey=ax1)
        semilogy(tt_relax*1.e3,T_aux_relax[:,0], label='Tx',nonposy='mask')
        semilogy(tt_relax*1.e3,T_aux_relax[:,1], label='Ty',nonposy='mask')
        semilogy(tt_relax*1.e3,T_aux_relax[:,2], label='Tz',nonposy='mask')
        semilogy(t_aux1*1.e3,y1_aux,'r')
        semilogy(t_aux2*1.e3,y1_aux,'r')
        ax2.set_yticks([1e-4,1e-2,1,1e2,1e4])
        ax2.grid()
        # annotate('laser off après injection', xy=(0.5,350), xycoords='data',
            # size=24, ha='left', va='top', color='xkcd:azul',
            # bbox=dict(boxstyle='round', fc='white',edgecolor='xkcd:azul'))
        legend(title='Laser OFF après injection',fontsize=18)    
        plt.setp(ax2.get_xticklabels(),visible=False)

        ax3 = subplot(313,sharex=ax1)
        plot(tt*1.e3,PM[:])
        plot(t_aux1*1.e3,y2_aux,'r')
        plot(t_aux2*1.e3,y2_aux,'r')
        ax3.grid()
        
        xlabel('time[ms]')
        ylabel('Counts')
                               
        ax1.set_xlim(xlim1)
        ax1.set_ylim(ylim1)
        ax3.set_ylim(ylim3)
        plt.tight_layout()
        
        fig.set_size_inches(11.69, 8.27)
        subplots_adjust(hspace=0.015)
        
        temps = [tt,t_aux1, t_aux2]
        fluo = [PM, PM_variation, SNR]
        temperature = [T_aux, T_variation, T_CM]
        
    return temps, temperature, fluo

def plot_T_and_PM_InitQ_Evol_AfterCool(data, flag_plot, fig_name, **kwargs):
    onlyfiles = sort([f for f in listdir(data) if isfile(join(data, f)) and "SimuType" in f and ".dat" in f])
    onlyfiles_Lan = sort([f for f in listdir(data) if isfile(join(data, f)) and "Harmo" in f])

    xlim1 = kwargs.get('xlim1', (-0.1, 6))
    ylim1 = kwargs.get('ylim1', (0.5 * 1e-3, 2e4))
    ylim2 = kwargs.get('ylim2', (-2, 50))

    N_ions, j_save, dt_j_save_next, eta, Temp, save_T = \
        load_Temp_init_bin_Lan(data + os.path.sep + onlyfiles_Lan[0].strip('.bin'), 0)
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(data + os.path.sep + onlyfiles[1].strip('.dat'))
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(data + os.path.sep + onlyfiles[0].strip('.dat'))

    tt0 = save_T[:, 0]
    T_CM0 = save_T[:, 1:4]
    T_aux0 = save_T[:, 4:]
    tt = concatenate((tt0, tt1, tt3))
    T_CM = concatenate((T_CM0, T_CM1, T_CM3))
    T_aux = concatenate((T_aux0, T_aux1, T_aux3))

    if flag_plot == 1:
        # fig_name = file_name[-9:]
        fig = figure(fig_name);
        clf()
        ax1 = subplot(211)
        semilogy(tt * 1.e3, T_aux[:, 0], label='Tx')
        semilogy(tt * 1.e3, T_aux[:, 1], label='Ty')
        semilogy(tt * 1.e3, T_aux[:, 2], label='Tz')
        ax1.grid()

        legend()
        # ~ xlabel('time[ms]')
        # ~ ylabel('T[K]')
        plt.setp(ax1.get_xticklabels(), visible=False)

        ax2 = subplot(212, sharex=ax1)
        plot(tt * 1.e3, T_aux[:, 0], label='Tx')
        plot(tt * 1.e3, T_aux[:, 1], label='Tx')
        plot(tt * 1.e3, T_aux[:, 2], label='Tx')
        ax2.grid()

        xlabel('time[ms]')
        ylabel('T [k]')

        ax1.set_xlim(xlim1)
        ax1.set_ylim(ylim1)
        ax2.set_ylim(ylim2)
        plt.tight_layout()
        subplots_adjust(hspace=0.015)

    return tt, T_aux, T_CM

def find_PM_variation_FinalT(file_dir2,file_name):
    i_aux = file_name.find('_N')+1
    file1 = 'SimuTypeQ_{}'.format(file_name[i_aux:].strip('50eV.dat'))    # for N=1024
    file3 = 'SimuType2_01_{}50eV'.format(file_name[i_aux:].strip('50eV.dat'))
    
    tt1, T_CM1, T_aux1, PM1 = load_T_and_PM_simu(file_dir2+'/Temp_'+file1)
    tt3, T_CM3, T_aux3, PM3 = load_T_and_PM_simu(file_dir2+'/Temp_'+file3)
    
    # 1ms is 2000 pts at 1pt/500ns
    
    aux = mean(PM1[-2000:])
    PM_variation = ( aux - mean(PM3[-2000:]) ) / aux
    T_variation  = ( mean(T_aux3[-2000:,0]) + mean(T_aux3[-2000:,1]) + mean(T_aux3[-2000:,2]) ) / 3
    
    SNR = np.abs( aux - mean(PM3[-2000:]) )/np.sqrt(aux) # Sig-to-Noi ratio considering Poisson noise
    
    return PM_variation, T_variation, SNR

# Look at the lost of energy:
def energy_lost(file_dir2,file_name):
    
    data = loadtxt('{}/{}.dat'.format(file_dir2,file_name))

    r_0 = data[0:3]; v_0 = data[3: 6]; # Position and velocity at the start
    r_1 = data[6:9]; v_1 = data[9:12]; # Position and velocity at the end
    t_c = data[12] ; # Time of crossing

    Ec0 = sum(v_0*v_0) ; Ec1 = sum(v_1*v_1)
    delta_Ec = (Ec1 - Ec0)
    return r_0, r_1, delta_Ec, delta_Ec / Ec1, t_c
    
def load_trj(str_load):
   
    fid = open(str_load+'.bin', 'rb')

    junk = fromfile(fid, int32,1)        # Read record start tag
    aux  = fromfile(fid, int32,1)
    junk = fromfile(fid, int32,1)        # Read record stop tag
    jmax = aux[0]

    junk = fromfile(fid, int32,1)        # Read record start tag
    aux  = fromfile(fid, int32,1)
    junk = fromfile(fid, int32,1)        # Read record stop tag
    N_ions1 = aux[0]

    junk = fromfile(fid, int32,1)        # Read record start tag
    aux  = fromfile(fid, int32,1)
    junk = fromfile(fid, int32,1)        # Read record stop tag
    N_ions2 = aux[0]
 
    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    t    = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    r_x  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    r_y  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    r_z  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    v_x  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    v_y  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    v_z  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    a_x  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    a_y  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    junk = fromfile(fid, int32  ,1   )        # Read record start tag
    a_z  = fromfile(fid, float64,jmax)
    junk = fromfile(fid, int32  ,1   )        # Read record stop tag

    fid.close

    return t,r_x, r_y, r_z, v_x, v_y, v_z, a_x, a_y, a_z

def load_cloud_trj(str_load):
    fid = open(str_load + '.bin', 'rb')

    junk0 = fromfile(fid, int32, 1)  # Read record start tag
    aux = fromfile(fid, int32, 1)
    junk1 = fromfile(fid, int32, 1)  # Read record stop tag
    jmax = aux[0]

    junk2 = fromfile(fid, int32, 1)  # Read record start tag
    aux = fromfile(fid, int32, 1)
    junk3 = fromfile(fid, int32, 1)  # Read record stop tag
    N_ions = aux[0]

    junk4 = fromfile(fid, int32, 1)  # Read record start tag
    r = fromfile(fid, single, int(junk4[0] / 4))  # jmax*N_ions*3
    junk5 = fromfile(fid, int32, 1)  # Read record stop tag

    fid.close

    return jmax, N_ions, r
    
##### FORTRAN ANALYSIS #####
###   Fonctions Adrien   ###

def load_gui(filter_nocomplete):

    # sélection des fichiers (Python 3)
    # SELECTIONNER UN FICHIER DE DONNÉES AU FOND DE L'ARBORESCENCE
    # Temp_SimuType0_N01024_Vrf0064_Udc0.5000D+00V_D1.0_S1.0RFG.dat

#     root = tk.Tk()
#     root.withdraw()

#     # Sélectionner un fichier du dossier pour récupérer le répertoire et les noms
#     file_path = filedialog.askopenfilename(initialdir = '/home/adrian/Documents/Programmes/Fortran_Jofre',
#                                            multiple=None)
    
    file_path = load_file_GUI('/home/adrian/Documents/Programmes/Fortran_Jofre')[0]
    
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
            slashcond = 6 # the slash number just before runs (after date)
        else:
            slashcond = 6 # the slash number just before runs (after date)
    else:
        slashcond = -2

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

def load_Temp_init_bin_Lan(str_load, flag_print):      

    fid = open(str_load+'.bin',  'rb')   

    a    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  int32, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    N_ions = aux[0]   

    b    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  int32, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    j_save = aux[0]   

    c    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  float64, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    dt_j_save_next = aux[0]   

    d    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  float64, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    eta = aux[0]   

    e    = fromfile(fid,  int32,   1   )        # Read record start tag   
    Temp = fromfile(fid,  float64, 1)   
    junk = fromfile(fid,  int32,   1   )        # Read record stop tag   

    f    = fromfile(fid,  int32,   1   )        # Read record start tag 
    save_T = fromfile(fid,  float64, 7*j_save)   
    junk = fromfile(fid,  int32,   1   )        # Read record stop tag

    fid.close   
    
    # print('len save_T', junk)   

    save_T = reshape(save_T, (j_save, 7), order='F')   
    
    if flag_print == 1:
        print(a, b, c, d, e, f)   
        print('N_ions', N_ions)   
        print('j_save', j_save)   
        print('dt_j_save_next', dt_j_save_next)   
        print('eta', eta)   
        print('Temp', Temp) 
        print(save_T[0])

    return N_ions, j_save, dt_j_save_next, eta, Temp, save_T 

def load_xva_init_bin_Lan(str_load, flag_print):      

    fid = open(str_load+'.bin',  'rb')   

    a    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  int32, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    N_ions = aux[0]   

    b    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  int32, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    j_save = aux[0]   

    c    = fromfile(fid,  int32, 1)        # Read record start tag   
    t_act= fromfile(fid,  float64, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   

    d    = fromfile(fid,  int32, 1)        # Read record start tag   
    dt   = fromfile(fid,  float64, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag
    
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

    return N_ions, j_save, t_act, dt, iRF, save_xva

#############################
### Single case functions ###

def file_name_retrieve(address,simu_type):
        # get only .dat files in each simulation directory
    onlyfiles = [f for f in listdir(address) if isfile(join(address, f)) and not "xva" in f and ".dat" in f]
    tmp_file = []
    # build path file
    if simu_type == 'Collision':
        for k in range(3):
            tmp_file.append('{}/{}'.format(address,sort(onlyfiles)[k].strip('.dat')))
    elif simu_type == 'RF_relax':
        for k in range(4):
            try:
                tmp_file.append('{}/{}'.format(address,sort(onlyfiles)[k].strip('.dat')))
            except(IndexError):
                tmp_file.append('No file for this SimuType')
                print('A file is missing, skiping the file number',k)
            
    return onlyfiles, tmp_file

# Return Temp, Fluo, for SimuType2 RF with laser
def load_T_PM_cloud_GMol(address,flag_plot):
    # retrieve files adresses
    onlyfiles, data_files = file_name_retrieve(address,simu_type='Collision')
    file0 = data_files[0];file2 = data_files[1];file4 = data_files[2]
    # load Ec variation for GMol
    _,_,deltaEc,deltaEcRel,t_c = energy_lost(address,'xva'+sort(onlyfiles)[2][4:].strip('.dat'))
    # load cloud size before injection
    r_LC_clip, dim_nu = cloud_size(address,onlyfiles)
    # load fluorescence and T var
    # Temperature during time
    temps, temperature, fluo = plot_T_and_PM_Init_Inje_Evol(address+'/',sort(onlyfiles)[0].strip('.dat')[4:],flag_plot,'T_PM',xlim1=(-1.5,55),ylim1=(5e-5,12e3))

    fluo_var = [fluo[1], temperature[1], fluo[2]]                        # temps = [tt,t_aux1, t_aux2]
    GMol_var = [deltaEc,deltaEcRel,t_c]                                  # fluo = [PM, PM_variation, SNR]
    cloud_atlas = [temps[0], temperature[0], temperature[2], fluo[0]]    # temperature = [T_aux, T_variation]
    t_aux = [temps[1],temps[2]]
    
    return fluo_var, GMol_var, cloud_atlas, t_aux, r_LC_clip
