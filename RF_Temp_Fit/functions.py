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
