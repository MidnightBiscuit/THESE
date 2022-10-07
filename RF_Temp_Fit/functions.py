import numpy as np
from matplotlib import pylab, mlab, pyplot
from pylab import *
from IPython.core.pylabtools import figsize, getfigs
plt = pyplot

import os
from os import listdir
from os.path import isfile, join
import sys

def beta_continue_alamano(a,q,beta_guess): #a l'air ok
    C1 = q**2/( (beta_guess+2)**2 - a - q**2/( (beta_guess+4)**2 - a - q**2/( (beta_guess+6)**2 - a - q**2/( (beta_guess+8)**2 - a - q**2/( (beta_guess+10)**2 - a ) ) ) ) )
    C2 = q**2/( (beta_guess-2)**2 - a - q**2/( (beta_guess-4)**2 - a - q**2/( (beta_guess-6)**2 - a - q**2/( (beta_guess-8)**2 - a - q**2/( (beta_guess-10)**2 - a  ) ) ) ) )
    return a + C1 + C2
    
def load_xv2(filename):
	fid = open(filename, 'rb')
	junk = fromfile(fid, int32,1)        # Read record start tag
	j_save_temp  = fromfile(fid, int32, 1);
	junk = fromfile(fid, int32,1)        # Read record stop tag

	junk = fromfile(fid, int32,1)        # Read record start tag
	n_ions  = fromfile(fid, int32, 1);
	junk = fromfile(fid, int32,1)        # Read record stop tag

	junk = fromfile(fid, int32,1)        # Read record start tag
	print(junk)
	r2_v2  = fromfile(fid, float64, junk[0]//8);
	junk = fromfile(fid, int32,1)        # Read record stop tag
	print(junk)

	fid.close

	r2_v2 = reshape(r2_v2, (6, junk[0]//8//6),order='F')
	return r2_v2
	
def load_xv2rlim(filename):
	fid = open(filename, 'rb')
	junk = fromfile(fid, int32,1)        # Read record start tag
	j_save_temp  = fromfile(fid, int32, 1);
	junk = fromfile(fid, int32,1)        # Read record stop tag

	junk = fromfile(fid, int32,1)        # Read record start tag
	n_ions  = fromfile(fid, int32, 1);
	junk = fromfile(fid, int32,1)        # Read record stop tag

	junk = fromfile(fid, int32,1)        # Read record start tag
	print(junk)
	r2_v2  = fromfile(fid, float64, junk[0]//8);
	junk = fromfile(fid, int32,1)        # Read record stop tag
	print(junk)

	fid.close

	r2_v2 = reshape(r2_v2, (9, junk[0]//8//9),order='F')
	return r2_v2

def load_x_3DHarmo(filename):
	fid = open(filename, 'rb')

	junk = fromfile(fid, int32,1)        # Read record start tag
	n_ions  = fromfile(fid, int32, 1);
	junk = fromfile(fid, int32,1)        # Read record stop tag

	junk = fromfile(fid, int32,1)        # Read record start 
	jj  = fromfile(fid, int32, 1);
	junk = fromfile(fid, int32,1)        # Read record stop tag

	junk = fromfile(fid, int32,1)        # Read record start tag
	dt  = fromfile(fid, float64, 1);
	junk = fromfile(fid, int32,1)        # Read record stop tag

	junk = fromfile(fid, int32,1)        # Read record start tag
	iRF  = fromfile(fid, int32, 3);
	print(iRF)
	junk = fromfile(fid, int32,1)        # Read record stop tag

	junk = fromfile(fid, int32,1)        # Read record start tag
	print(junk)
	xva_3DHarmo  = fromfile(fid, float64, junk[0]//8);
	junk = fromfile(fid, int32,1)        # Read record stop tag
	print(junk)

	fid.close

	xva_3DHarmo = reshape(xva_3DHarmo, (12, junk[0]//8//12),order='F')
	return xva_3DHarmo
	
def load_x_afterstart(filename):
	fid = open(filename, 'rb')

	junk = fromfile(fid, int32,1)        # Read record start tag
	xva_trj1  = fromfile(fid, float64, junk[0]//8);
	junk = fromfile(fid, int32,1)        # Read record stop tag

	fid.close

	xva_trj1 = reshape(xva_trj1, (6, junk[0]//8//6),order='F')
	return xva_trj1

def load_xyz_init_bin_DP(str_load):
    fid = open(str_load+'.bin',  'rb')   
    a    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  int32, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    N_ions = aux[0]   

    b    = fromfile(fid,  int32, 1)        # Read record start tag   
    aux  = fromfile(fid,  int32, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag   
    j_save = aux[0]     

    d    = fromfile(fid,  int32, 1)        # Read record start tag   
    dt   = fromfile(fid,  float64, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag
    
    e    = fromfile(fid,  int32,   1   )        # Read record start tag   
    iRF  = fromfile(fid,  int32, 3)   
    junk = fromfile(fid,  int32,   1   )        # Read record stop tag   

    f    = fromfile(fid,  int32,   1   )        # Read record start tag 
    # ~ print(f)
    save_xva = fromfile(fid,  float64, 12*N_ions)   
    junk = fromfile(fid,  int32,   1   )        # Read record stop tag
    # ~ print(junk)

    fid.close
    save_xva = reshape(save_xva, (12,N_ions), order='F')

    return save_xva

def plot_XYZ(file_name,fig_name='2',fig_title='XYZ'):
    save_xva = load_xyz_init_bin_DP(file_name)
    r_LC = save_xva[0:3,:] 
    v_LC = save_xva[3:6,:]
    a_LC = save_xva[6:9,:]
    v_LC_avg = save_xva[9:,:]
    # ~ print(r_LC)
    figure(fig_name); clf()
    title(fig_title)
    subplot(211)
    plot(r_LC [0,:]*1e6,r_LC [1,:]*1e6,'8',color='xkcd:purplish blue')
    xlabel('x [µm]')
    ylabel('y [µm]')
    grid()

    # subplot(212,aspect=1.0)
    subplot(212,aspect=1.0)
    plot(r_LC [2,:]*1e6,r_LC [1,:]*1e6,'8',color='xkcd:purplish blue')
    xlabel('z [µm]')
    ylabel('y [µm]')
    grid()
    
    tight_layout()

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

    d    = fromfile(fid,  int32, 1)        # Read record start tag   
    dt   = fromfile(fid,  float64, 1)   
    junk = fromfile(fid,  int32, 1)        # Read record stop tag
    
    e    = fromfile(fid,  int32,   1   )        # Read record start tag   
    iRF  = fromfile(fid,  int32, 3)   
    junk = fromfile(fid,  int32,   1   )        # Read record stop tag   

    f    = fromfile(fid,  int32,   1   )        # Read record start tag 
    # ~ print(f)
    save_xva = fromfile(fid,  float64, 12*N_ions)   
    junk = fromfile(fid,  int32,   1   )        # Read record stop tag
    # ~ print(junk)

    fid.close
    save_xva = reshape(save_xva, (12,N_ions), order='F')

    return N_ions, j_save, dt, iRF, save_xva
