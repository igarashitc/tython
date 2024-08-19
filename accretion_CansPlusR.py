import numpy as np
from scipy import constants as scp
import matplotlib.pyplot as plt

from dac_read_mpiio import dac_read
from util_CansPlusR import util_cyl
import plot_CansPlusR as plcans

class accretion(util_cyl):

  #=========================
  # pysical constatns in cgs
  #
  # speed of light
  cc = scp.c*1e2
  #boltzmann constant
  kb = scp.k*1e7
  #stefan-boltzmann constant
  sb = scp.sigma*1e3
  #radiation constant
  ar = 4*sb/cc
  #proton mass 
  mp = scp.m_p*1e3
  #electron mass
  me = scp.m_e*1e3
  #gravitational constant
  gc = scp.G*1e3
  #solar mass
  Ms = 1.98841586057e33
  #=========================

  def __init__(self, num, dir_name, density=1e0, black_hole=1e1, mean_molecular_weight=0.5, *args):

    #===========#
    #parameters #
    #===========#
    #normalized density
    self.de = density
    #black hole mass
    self.bh = black_hole
    #schwarzchild radius
    self.rs = 2e0*self.gc*self.bh*self.Ms/(self.cc*self.cc)
    #mean molecular_weight
    self.mmw = mean_molecular_weight
    #===============#
    
    #===============#
    #physical value #
    #===============#
    #electron scattering opacity
    self.kes = 0.4
    #temperature
    self.te0 = self.mmw*self.mp*self.cc*self.cc/self.kb
    #Eddington luminosity
    self.Ledd = 4*np.pi*self.cc*self.gc*self.bh*self.Ms/self.kes
    #Eddington accretion rate
    self.Mded = self.Ledd/(self.cc*self.cc)
    #===============#

    # data number
    self.num = int(num)

    # data input
    data_dir = dir_name+str(self.num).zfill(4)+"_"

    print("input dir and data number:",data_dir)
    self.val = {}

    tmp,r,p,z = dac_read(data_dir+"ro"+".dac")
    self.val["mass_dens"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"pr"+".dac")
    self.val["gas_press"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"vx"+".dac")
    self.val["radial_vel"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"vy"+".dac")
    self.val["azimuthal_vel"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"vz"+".dac")
    self.val["vertical_vel"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"bx"+".dac")
    self.val["radial_mag"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"by"+".dac")
    self.val["azimuthal_mag"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"bz"+".dac")
    self.val["vertical_mag"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"er"+".dac")
    self.val["rad_ene"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"fx"+".dac")
    self.val["radial_rad_flux"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"fy"+".dac")
    self.val["azimuthal_rad_flux"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"fz"+".dac")
    self.val["vertical_rad_flux"] = tmp[0,:,:,:].astype(np.float16)

    print("Successfully loaded.",self.val.keys())
    #coordinate in cylindrical
    if (len(args) == 0):
      super().__init__(r,p,z)
    else:
      super().__init__(r,p,z,args)

    self.fig = plt.figure()

  #======================
  # reload data
  #======================
  def load_val(self, num):
    # data number
    self.num = int(num)

    # data input
    data_dir = dir_name+str(self.num).zfill(4)+"_"

    print("input dir and data number:",data_dir)
    self.val = {}

    tmp,r,p,z = dac_read(data_dir+"ro"+".dac")
    self.val["mass_dens"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"pr"+".dac")
    self.val["gas_press"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"vx"+".dac")
    self.val["radial_vel"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"vy"+".dac")
    self.val["azimuthal_vel"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"vz"+".dac")
    self.val["vertical_vel"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"bx"+".dac")
    self.val["radial_mag"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"by"+".dac")
    self.val["azimuthal_mag"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"bz"+".dac")
    self.val["vertical_mag"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"er"+".dac")
    self.val["rad_ene"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"fx"+".dac")
    self.val["radial_rad_flux"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"fy"+".dac")
    self.val["azimuthal_rad_flux"] = tmp[0,:,:,:].astype(np.float16)
    tmp,r,p,z = dac_read(data_dir+"fz"+".dac")
    self.val["vertical_rad_flux"] = tmp[0,:,:,:].astype(np.float16)

    print("Successfully loaded.",val.keys())

  #======================
  # reset figure
  #======================
  def set_fig(self):

    del self.fig
    self.fig = plt.figure()

#==================#
# calcurate values #
#==================#

  #===============================
  # calcurate optical 
  #  depth for electron scattering.
  # calcurate from over 
  #  the disk to equatorial plane.
  #===============================
  def electron_optical_depth_z(self):

    tau = np.zeros([self.nz,self.ny,self.nr], np.float16)

    z0 = np.argmin(abs(self.z))

    ro = self.val["mass_dens"]

    for i in range(self.nr):
      for j in range(self.ny):
        for k in range(1, z0):
          tau[k,j,i] = tau[k-1,j,i] + ro[k,j,i]*self.dz[k]
        for k in range(self.nz-1, z0+1):
          tau[k,j,i] = tau[k,j,i] + ro[k,j,i]*self.dz[k]
    
    self.val["electron_optical_depth"] = tau

    return tau

  #===============================
  # calcurate surface density
  #===============================
  def vert_int_value(self, *args):

    sig = super().vertical_integrate(self.val["mass_dens"], args)

    return sig
  
  #===============================
  # calcurate vertically 
  #  integrated values.
  #===============================
  def surface_density(self, data_name, *args):

    if( data_name in self.val ):
      print("data name missmach.")
      print(self.val.keys())
      exit

    sig = super().vertical_integrate(self.val[data_name], args)

    return sig

  #===============================
  # calcurate accretion rate.
  #===============================
  def accretion_rate(self, *args):

    mass_flux = self.val["mass_dens"]*self.val["radial_vel"]
    mdot_net =-super().surface_integrate(mass_flux, args)
    mdot_out = super().surface_integrate(self.val["mass_dens"], args)
    mdot_in  =-super().surface_integrate(self.val["mass_dens"], args)

    return mdot_net, mdot_in, mdot_out

  #===============================
  # calcurate luminosity.
  #===============================
  def luminosity(self, *args):

    lu = np.ndarray(2)

    print("Luminosity ::", lu, "L_{Edd}")

    return lu

#================#
# plot procedure #
#================#
  
  #===============================
  # p0lar pl0t
  #===============================
  def plot_polar(self, data_name, scale=0, cmap="coolwarm", \
                  x1=-80, x2=80, y1=-80, y2=80):

    import matplotlib.colors as colors
    
    #====================
    # data setting
    data = self.val[data_name]
    tmp = self.p.copy()
    tmp[1:] = tmp[1:]+self.dp
    Rad,Phi = np.meshgrid(self.r-0.5*self.dr,tmp)
    X,Y = Rad*np.cos(Phi), Rad*np.sin(Phi)
    #====================
    
    #====================
    # linear or log scale
    if (scale == 0):
      norm = colors.Normalize(vmin=data.min(),vmax=data.max())
    else: 
      if (scale == 1):
        norm = colors.LogNorm(vmin=data.min(),vmax=data.max())
    #====================

    #====================
    # plot
    fig,ax = plt.add_subplots() 
    pc = ax.pcolormesh(X,Y,data[z0,:,:],norm=norm,cmap=cmap)
    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)
    ax.set_aspect("equal")
    #====================

    return fig,ax,pc
  
  #===============================
  # c0l0rbar
  #===============================
  def colorbar_single(self, pc, ax):

    cb = fig.colorbar(pc, ax=ax)

    return cb

  #===============================
  # pl0t slice in rz-plane.
  #===============================
  def plot_rzslice_single(self, data_name, scale=0, avr=1, \
                  r1=0, r2=80, z1=-20, z2=60):
    
    import matplotlib.colors as colors

    #====================
    # data setting
    if(avr==0):
      cdata = self.val[data_name]
    else:
      if (avr==1):
        cdata = average(self.val[data_name],axis=1)
    #====================
    
    #====================
    # linear or log scale
    if (scale == 0):
      norm = colors.Normalize(vmin=cdata.min(),vmax=cdata.max())
    else: 
      if (scale == 1):
        norm = colors.LogNorm(vmin=cdata.min(),vmax=cdata.max())
    #====================

    #====================
    # data setting
    fig,ax = plt.subplots()
    pc = ax.pcolormesh(self.r,self.z,cdata,norm=norm,cmap=cmap)
    ax.set_xlim(r1,r2)
    ax.set_ylim(z1,z2)
    ax.gca().set_aspect("equal")
    #====================

    return fig,ax,pc

  #===============================
  # pl0t slice in rz-plane.
  #=============================== 
  def plot_rzslice_multi(self, data_name_l, data_name_r,\
                        scale_l=0, scale_r=1, avr_l=0, avr_r=0,\
                        cmap_l="coolwarm", cmap_r="nipy_spectral",\
                        r1=-80, r2=80, z1=-40, z2=120):

    data_l = val["data_name_l"]
    data_r = val["data_name_r"]

    ax,pcl,pcr = plcans.plot_rzslice_multi(\
                          self.r, self.z,data_l, data_r, self.fig,\
                          scale_l=scale_l, scale_r=scale_r,\
                          avr_l=avr_l, avr_r=0,\
                          cmap_l=cmap_l, cmap_r=cmap_r,\
                          r1=r1,r2=r2,z1=z1,z2)

    return fig,pc_l,pc_r
  
