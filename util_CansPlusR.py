import numpy as np

class util_cyl:

  def __init__(self, r, p, z, *args):
    #coordinate in cylindrical
    self.r = r.astype(np.float16)
    self.p = p.astype(np.float16)
    self.z = z.astype(np.float16)

    #number of grid point
    self.nr = r.shape[0]
    self.ny = p.shape[0]
    self.nz = z.shape[0]

    #grid spacing
    if (len(args) == 0):
      self.dr = np.ndarray(self.nr, np.float16)
      for i in range(1,self.nr-1):
        self.dr[i] = (r[i+1]-r[i-1])*0.5

      self.dp    = np.ndarray(self.ny, np.float16)
      if (self.ny == 1):
        self.dp[0] = 2.0*np.pi
      else:
        self.dp[:] = self.p[1]-self.p[0]

      self.dz = np.ndarray(self.nz, np.float16)
      for k in range(1,self.nz-1):
        self.dz[k] = (z[k+1]-z[k-1])*0.5

    #given grid spacing
    if (len(args) != 0):
      self.dr = args[0].astype(np.float16)
      self.dp = args[1].astype(np.float16)
      self.dz = args[2].astype(np.float16)

  #integrate over vertical direction
  def integrate_dz(self, *args):
   
    val = args[0]

    srf = find_surface(self, args[1:])

    sig = np.zeros([self.ny, self.nr])

    for i in range(self.nr):
      for j in range(self.ny):
        for k in range(srf[0,j,i],srf[1,j,i]):
          sig[j,i] = sig[j,i] + val[k,j,i]*self.dz

    return sig

  def integrate_rdphidz(self, *args):

    val = args[0]

    srf = find_surface(self, args)

    sig = np.zeros([self.ny,self.nr])

    for i in range(self.nr):
      for j in range(self.ny):
        for k in range(srf[0,j,i],srf[1,j,i]):
          sig[i] = sig[i] + val[k,j,i]*self.r[i]*self.dp[j]*self.dz
   
    return sig
  
  def integrate_rdphidr(self, *args):

    val = args[0]

    srf = find_surface(self, args)

    sig = np.zeros([2,self.ny,self.nr])

    for j in range(self.ny):
      for i in range(self.nr):
        #upper side
        z1 = srf[0,i]
        sig[0,j,i] = sig[0,j,i] + val[z1,j,i]*self.r[i]*self.dp[j]*self.dz
        #lower side
        z2 = srf[1,i]
        sig[1,j,i] = sig[1,j,i] + val[z2,j,i]*self.r[i]*self.dp[j]*self.dz
   
    return sig

  def integrate_volume(self,*args,ri=0,rf=200,z1=-10,z2=10):

    val = args[0]

    z1p = np.argmin(np.abs(self.z+z1))
    z1m = np.argmin(np.abs(self.z-z2))
    z0  = np.argmin(np.abs(self.z))

    sf = find_surface(depth, z, zpt=zpt)

    r1 = np.argmin(np.abs(self.r-ri))
    r2 = np.argmin(np.abs(self.r-rf))

    vol = 0.0
    for i in range(r1,r2):
      for j in range(self.ny):
        for k in range(sf[0,j,i],sf[1,j,i]):
          vol = vol + val[k,j,i]*self.r[i]*self.dp[j]*self.dr[i]*self.dz[k]

    return vol
      
  #find optical depth = zpt
  def find_surface(self, *args):
    
    surf = np.zeros([2,self.ny,self.nr])

    z0 = np.argmin(np.abs(self.z))

    if (len(args) == 0):
      zpt = 40.0
      zu  = np.argmin(np.abs(self.z-zpt))
      zd  = np.argmin(np.abs(self.z+zpt))
      srf[0,:,:] = dz
      srf[1,:,:] = du
      print("integrate between +-", zpt)

    elif (type(args[0]) == np.float):
      zpt = args[0]
      zu  = np.argmin(np.abs(self.z-zpt))
      zd  = np.argmin(np.abs(self.z+zpt))
      srf[0,:,:] = zd
      srf[1,:,:] = zu
      print("integrate between +-", zpt)

    elif (type(args[0]) == np.ndarray):
      if (len(args) < 2 and type(args[1]) != np.float):
        print("please specify the point.")

      for j in range(self.ny):
        for i in range(self.nr):
          tmp = tau[z0,j,i]-zpt
          surf[0,j,i] = np.argmin(np.abs(tau[:z0,j,i]-zpt))
          surf[0,j,i] = surf[0,j,i] + \
                        (z0+1 - surf[0,j,i])*min(tmp, 0)/tmp

          tmp = tau[z0+1,j,i]-zpt
          surf[1,j,i] = max( \
                        np.argmin(np.abs(tau[z0+1:,j,i]-zpt))+z0+1,\
                        z0*max(tmp, 0)/tmp)
          
          surf = np.asarray(surf,int)

    return surf

#find optical depth = zpt
def surface_opt_av_azim(self, tau, zpt=10):

    surf = np.zeros([2,self.ny,self.nr])

    z0 = np.argmin(abs(self.z))
    mtau = np.average(tau,axis=1)

    for i in range(self.nr):
        tmp = mtau[z0,i]-zpt
        surf[0,:,i] = np.argmin(abs(mtau[:z0,i]-zpt)) 
        surf[0,:,i] = surf[0,:,i] + (z0+1 - surf[0,:,i])*min(tmp, 0)/tmp
        tmp = mtau[z0+1,i]-zpt
        surf[1,i] = max(np.argmin(abs(mtau[z0+1:,i]-zpt))+z0+1, z0*max(tmp, 0)/tmp)
    surf = np.asarray(surf,int)

    return surf

#class util_sphr():
  


