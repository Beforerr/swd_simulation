import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import pywt
from matplotlib import ticker
from scipy.fft import fft, fft2, fftfreq, fftshift
from scipy.signal import savgol_filter

datadir = "/glade/derecho/scratch/xinan/low2high/hvpic/1/data/"
nx = 2048#96
ny = 256
nt = 81
slice = int(input('type time stamps: '))

pi = np.pi

xv = np.linspace(0,120,num=nx)
yv = np.linspace(0,15,num=ny)
tv = np.linspace(0,800,num=nt)
if (nx>1): dx = xv[1]-xv[0]
if (ny>1): dy = yv[1]-yv[0]
if (nt>1): dt = tv[1]-tv[0]

########## wavelet denoising function
def wclean(arr,wavn,alpha):
	cs = pywt.wavedecn(arr, wavn, mode='symmetric', level=None,axes=None)
	levs = len(cs)
	coef = np.concatenate(cs[0])
	for x in range(1,levs):
		for n in cs[x]:
			coefn = np.concatenate(cs[x][n])
			coef = np.concatenate((coef,coefn))	
	thr2 = alpha*np.sqrt(np.var(coef)*np.log(len(coef)))
	thr = alpha*2*thr2
	while (thr/thr2 > 1.05):
		thr = thr2;
		thr2 = alpha*np.sqrt(np.var(coef[np.abs(coef)<thr])*np.log(len(coef)))
	for x in range(1,levs):
		for n in cs[x]:
			inds = np.abs(cs[x][n]) < thr
			cs[x][n][inds] = 0
	return pywt.waverecn(cs, wavn, mode='symmetric', axes=None);
########## end wavelet denoising


######### loadSlice function
def loadSlice(dir,q,sl,nx,ny):
	fstr = dir + q + ".gda"
	fd = open(fstr,"rb")
	fd.seek(4*sl*nx*ny,1)
	arr = np.fromfile(fd,dtype=np.float32,count=nx*ny)
	fd.close
	arr = np.reshape(arr,( ny, nx))
	arr = np.transpose(arr)
	return arr
######### end loadSlice

#yv,xv = np.meshgrid(np.linspace(0,7.5*pi,num=ny),
#                np.linspace(0,5*pi,num=nx))


Q = {}


# for slice in range(0,1):
# slice = 400
qs = ["ni","uix", "uiy", "uiz", "bz", "by", "bx",\
    "ex", "ey", "ez"]
for q in qs:
	tmp = loadSlice(datadir,q,slice,nx,ny)
	Q[q] = tmp
		
# Qplts = ['By', 'Bz', 'Uiy', 'Uiz', 'Uix', 'ni']
Qplts = ['ey', 'ex', 'uix', 'ni', 'bx', 'by']

cmap = plt.get_cmap("RdBu_r")

for j, qplt in enumerate(Qplts):
  fig, ax = plt.subplots(1,1, figsize=(6.4, 4.8))
  
  vmax = np.max(np.abs(Q[qplt]))
  if qplt == 'ni' or qplt == 'bx':
    vmin = 2 - vmax
  else:
    vmin = -vmax
  # vmin = np.min(Q[qplt])
  im = ax.pcolormesh(yv,xv,Q[qplt], vmin=vmin, vmax=vmax, cmap=cmap)

  ax.set_aspect('equal')
  
  # colorbar
  l = ticker.AutoLocator()
  l.create_dummy_axis()
  ticks = l.tick_values(vmin, vmax)
  cb = fig.colorbar(im, ax=ax, ticks=ticks, orientation='vertical')
  cb.ax.set_ylabel(qplt)
        
  ax.set_ylabel(r'$x/d_i$')
  ax.set_xlabel(r'$y/d_i$')
  # plt.show()
  plt.tight_layout()
  plt.savefig(datadir+qplt+'-'+str(slice)+'.png')
  plt.close()


# cmap = plt.get_cmap("nipy_spectral")
cmap = plt.get_cmap("BuPu")

for j, qplt in enumerate(Qplts):
  fig, ax = plt.subplots(1,1, figsize=(6.4, 4.8))
  # qf = fft(Q[qplt])[:nx//2,:ny//2]
  qf = fft2(Q[qplt])
  qf = fftshift(qf)
  qpow = (np.abs(qf))**2
  wavenumx = 2.0* np.pi * fftfreq(np.shape(xv)[0], d=dx)
  wavenumy = 2.0* np.pi * fftfreq(np.shape(yv)[0], d=dy)
  wavenumx = fftshift(wavenumx)
  wavenumy = fftshift(wavenumy)
  vmax = np.nanmax(np.log10(qpow))
  vmin = vmax - 7
  im = ax.pcolormesh(wavenumy,wavenumx,np.log10(qpow), vmin=vmin, vmax=vmax, cmap=cmap)
  
  # ax.set_xscale('log')
  ax.set_xscale('symlog', linthresh=1)
  ax.set_yscale('log')
  if qplt == 'ex' or qplt == 'ni' or qplt == 'uix':
    ax.set_ylim([0.01, 50])
    ax.set_xlim([-30, 30])
  else:
    ax.set_ylim([0.01, 50])
    ax.set_xlim([-30, 30])
  # ax.set_aspect('equal')
  
  # colorbar
  l = ticker.AutoLocator()
  l.create_dummy_axis()
  ticks = l.tick_values(vmin, vmax)
  cb = fig.colorbar(im, ax=ax, ticks=ticks, orientation='vertical')
  cb.ax.set_ylabel(r'$\log_{10}$'+qplt)
  ax.set_ylabel(r'$k_x d_i$')
  ax.set_xlabel(r'$k_y d_i$')
  plt.tight_layout()
  plt.savefig(datadir+qplt+'-spec'+'-'+str(slice)+'.png')
  plt.close()
