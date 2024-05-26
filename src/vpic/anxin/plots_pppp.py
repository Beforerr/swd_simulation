import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.signal import savgol_filter

datadir = "/glade/derecho/scratch/xinan/low2high/hvpic/0/data/"
ny = 2048
nu = 200
ymax = 120.0
umax = 2.0 
betai = 0.1
vti = np.sqrt(betai * 0.5)
dt = 0.004
# umax = 1.2 # run 12
# umax = 1.5 # run 11

# prefix = "xux.240000"
ps = 'xux'
sl = input("type prefix of particle data file: ")
prefix = ps + '.' + sl
data = np.fromfile(datadir+prefix+".bin", dtype=np.float32)
data = np.reshape(data, (nu,ny), order='C') 

# initial distribution function
xarr = np.linspace(0, ymax, ny)
uarr = np.linspace(-umax, umax, nu)
fv = np.exp(-0.5*(uarr/vti)**2)
fv = fv / np.sum(fv) * (np.sum(data)/ny)
f0 = np.transpose(np.tile(fv, (ny, 1)))

# plot
plt.style.use('~/python_lib/plt_style.txt')

fig, ax = plt.subplots(1, 1, figsize=(9.6, 3.2))

vmax = np.max(np.log10(data))
vmin = vmax - 3
cmap = plt.get_cmap("nipy_spectral")
# cmap = plt.get_cmap("BuPu")
# cmap = plt.get_cmap("Reds")
# cmap = plt.get_cmap("viridis")

im = ax.pcolormesh(xarr, uarr, np.log10(data), vmin=vmin, vmax=vmax, cmap=cmap)

# data_sm = savgol_filter(data,     63, 3, mode='nearest', axis=0)
# data_sm = savgol_filter(data_sm, 255, 3, mode='nearest', axis=1)
# levels = np.linspace(vmin, vmax, 8)
# cs = ax.contour(xarr, uarr, np.log10(data_sm), levels=levels,\
#     colors='k', linestyles='solid')

ax.set_xlim([xarr[0], xarr[-1]])
ax.set_ylim([uarr[0], uarr[-1]])
ax.set_xlabel(r'$x / d_i$')
ax.set_ylabel(r'$v_x / v_\mathrm{A}$')
# ax.set_ylabel(r'$v_\perp / v_\mathrm{A}$')

ax.text(0.5, 1.1, r"$t = {:.1f} \omega_{{ci}}^{{-1}}$".format(int(sl)*dt),
        horizontalalignment='center', verticalalignment='center',
        transform=ax.transAxes)

# colorbar
l = ticker.AutoLocator()
l.create_dummy_axis()
ticks = l.tick_values(vmin, vmax)
cb = fig.colorbar(im, ax=ax, ticks=ticks, orientation='vertical')
cb.ax.set_ylabel(r'$\log_{10}f$')

plt.tight_layout()
# plt.show()
plt.savefig(datadir+ps+'.{:06d}.png'.format(int(sl)))
plt.close()


# delta f
fig, ax = plt.subplots(1, 1, figsize=(9.6, 3.2))

deltaf = data - f0
vmax = np.max(np.abs(deltaf))
vmin = -vmax
# cmap = plt.get_cmap("nipy_spectral")
# cmap = plt.get_cmap("BuPu")
# cmap = plt.get_cmap("Reds")
# cmap = plt.get_cmap("viridis")
cmap = plt.get_cmap("RdBu_r")

im = ax.pcolormesh(xarr, uarr, deltaf, vmin=vmin, vmax=vmax, cmap=cmap)

ax.set_xlim([xarr[0], xarr[-1]])
ax.set_ylim([uarr[0], uarr[-1]])
ax.set_xlabel(r'$x / d_i$')
ax.set_ylabel(r'$v_x / v_\mathrm{A}$')
# ax.set_ylabel(r'$v_\perp / v_\mathrm{A}$')

# colorbar
l = ticker.AutoLocator()
l.create_dummy_axis()
ticks = l.tick_values(vmin, vmax)
cb = fig.colorbar(im, ax=ax, ticks=ticks, orientation='vertical')
cb.ax.set_ylabel(r'$\delta f$')

plt.tight_layout()
# plt.show()
plt.savefig(datadir+prefix+'.deltaf.png')
plt.close()

