
import numpy as np
import scipy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from util_CansPlusR import util_cyl

#==================
# polar plot
#==================
def plot_polar(r, p, dr, dz, val, *args,\
              scale="lin", cmap="coolwarm",\
              x1=-80, x2=80, y1=-80, y2=80,):

  #====================
  # data setting
  tmp = self.p.copy()
  tmp[1:] = tmp[1:]+self.dp
  Rad,Phi = np.meshgrid(self.r-0.5*self.dr,tmp)
  X,Y = Rad*np.cos(Phi), Rad*np.sin(Phi)
  #====================
  
  #====================
  # linear or log scale
  if (scale == "lin"):
    norm = colors.Normalize(vmin=val.min(),vmax=val.max())
  else:
    if (scale == "log"):
      norm = colors.LogNorm(vmin=val.min(),vmax=val.max())
  #====================

  #====================
  # plot
  if (len(args) == 0):
    ax = fig.add_subplot()
  else:
    ax = args[0]

  pc = ax.pcolormesh(X,Y,self.val,cmap=cmap,vmax=vmax,vmin=vmin,shading="gourand")
  pc = ax.set_aspect("equal")
  #====================

  return ax,pc
 
#====================
# plot poloidal plane
#====================
def plot_rzslice_multi(r, z, val_l, val_r, fig, *args,\
                  scale_l="lin", scale_r="log", avr_l=0, avr_r=0,\
                  cmap_l="coolwarm", cmap_r="nipy_spectral",\
                  r1=80., z1=-20., z2=60.):
    
  from mpl_toolkits.axes_grid1.inset_locator import inset_axes

  #=====================================
  # data setting
  # left
  if ( avr_l == 0 ):
    cval_l = val_l
  else:
    if ( avr_l == 1 ):
      cval_l = np.average(val_l,axis=1)
  # right
  if ( avr_r == 0 ):
    cval_r = val_r
  else:
    if ( avr_r == 1 ):
      cval_r = np.average(val_r,axis=1)

  tmp = np.fliplr(cval_l).copy()
  rf  = -np.flip(r)
  #=====================================
  
  #=====================================
  # linear or log scale
  # left
  if ( scale_l == "lin" ):
    norm_l = colors.Normalize(vmin=cval_l.min(),vmax=cval_l.max())
  else:
    if( scale_l == "log" ):
      norm_l = colors.LogNorm(vmin=cval_l.min(),vmax=cval_l.max())
  # right
  if ( scale_r == "lin" ):
    norm_r = colors.Normalize(vmin=cval_r.min(),vmax=cval_r.max())
  else:
    if ( scale_r == "log" ):
      norm_r = colors.LogNorm(vmin=cval_r.min(),vmax=cval_r.max())
  #=====================================

  #=====================================
  # plot 
  if ( len(args) == 0 ):
    ax = fig.add_subplot()
  else:
    ax = args[0]

  plt.rcParams["font.size"] = 14
  pc_l = ax.pcolormesh(rf, z, tmp   , cmap=cmap_l, norm=norm_l)
  pc_r = ax.pcolormesh(r , z, cval_r, cmap=cmap_r, norm=norm_r)
  ax.set_xlim(-r1,r1)
  ax.set_ylim( z1,z2)
  ax.set_aspect("equal")
  #=====================================

  #=====================================
  # colorbar
  # left
  axins_l = inset_axes(ax, width="45%", height="5%", loc="lower left",\
                      bbox_transform=ax.transAxes,\
                      bbox_to_anchor=(0,1.025,1,1))
  cb_l    = fig.colorbar(pc_l, cax=axins_l, orientation="horizontal")
  axins_l.xaxis.set_ticks_position("top")
  axins_l.xaxis.set_label_position("top")
  # right
  axins_r = inset_axes(ax, width="45%", height="5%", loc="lower right",\
                      bbox_transform=ax.transAxes,\
                      bbox_to_anchor=(0,1.025,1,1))
  cb_r    = fig.colorbar(pc_r, cax=axins_r, orientation="horizontal")
  axins_r.xaxis.set_ticks_position("top")
  axins_r.xaxis.set_label_position("top")
  #=====================================

  return ax, pc_l, pc_r

def set_crange(vmax, vmin, pc):

  pc.set_clim(vmax, vmin)


#=====================================
# plot psd at each radius.
#=====================================
def psd_scipy_2d(val, time, dt):
  import scipy.fftpack as sf
  
  nx   = val.shape[1]
  
  N    = len(time)
  posN = len(time[0:int(N/2)])
  
  psd  = np.zeros([posN, nx])
  
  freq = sf.fftfreq(N, dt)[0:int(N/2)]
  df   = freq[1]-freq[0]
  
  for i in range(nx):
    fft      = sf.fft(val[:,i], N)[0:int(N/2)]
    psd[:,i] = np.abs(fft)/np.sqrt(df)/(N/np.sqrt(2))
  
  return psd, freq
    
