import numpy as np

class accretion:

  def __init__(self, val, r, p, z):
    #value
    self.val = val

    #coordinate in cylindrical
    self.r = r
    self.p = p
    self.z = z

    self.nr = r.shape[0]
    self.ny = p.shape[0]
    self.nz = z.shape[0]

  def vertical_integrate(self, *args):
    
    sig = np.zeros([   self.nr])
    srf = np.zeros([2, self.nr])

    if (len(args) == 0):
      zpt = 40.0
      zu  = np.argmin(np.abs(self.z-zpt))
      zd  = np.argmin(np.abs(self.z+zpt))
      srf[0,:] = zd
      srf[1,:] = zu

    elif (type(args[0]) == np.float):
      zpt = args[0]
      zu  = np.argmin(np.abs(self.z-zpt))
      zd  = np.argmin(np.abs(self.z+zpt))
      srf[0,:] = zd
      srf[1,:] = zu

    elif (type(args[0]) == np.ndarray):
      print("on going")

    else:
      print("type missmach")

    srf = np.asarray(srf, np.int)
    print(srf[0,100],srf[1,100])

    for i in range(self.nr):
      for j in range(self.ny):
        for k in range(srf[0,i],srf[1,i]):
          dz = np.abs(self.z[k+1]-self.z[k-1])*0.5
          sig[i] = sig[i] + self.val[k,j,i]*dz

    return sig

      

def calc(ro,vx,x,y,z,density=1,length=1):
    from tython import radiation as rd

    nt = ro.shape[0]
    nx = ro.shape[3]

    mdot = np.ndarray(nt*nx).reshape(nt,nx)
    sig  = np.ndarray(nt*nx).reshape(nt,nx)
    
    for it in range(nt):
        tde = rd.thomson_depth(ro[it,:,:,:],z,density=density,length=length)
        #mdot[it,:] = mdot_net(ro[it,:,:,:],vx[it,:,:,:],tde,x,y,z)
        sig[it,:] = surface_density(ro[it,:,:,:],tde,x,y,z)
        print(it)

    return mdot,sig

def mdot(ro,vx,depth,x,y,z,zpt=1):
    nx = ro.shape[2]
    ny = ro.shape[1]
    nz = ro.shape[0]

    mdot_net = np.zeros(nx)
    mdot_out = np.zeros(nx)
    mdot_in  = np.zeros(nx)

    dy = abs(y[1]-y[0])

    sf = surface_opt_av_azim(depth, z, zpt=zpt)
    ## if the tau=zpi dose not exist
    z0 = np.argmin(abs(z))
    #print(30)
    for j in range(ny):
        for i in range(nx):
            if(abs(z[sf[0,i]]) < 0.5):
                sf[0,i] = np.argmin(abs(z+x[i]))
            if(abs(z[sf[1,i]]) < 0.5):
                sf[1,i] = np.argmin(abs(z-x[i]))

    for i in range(nx-1):
        for j in range(ny):
            for k in range(sf[0,i],sf[1,i]):
                dz = (z[k+1]-z[k-1])*0.5
                mdot_net[i] = mdot_net[i] - ro[k,j,i]*vx[k,j,i]*x[i]*dy*dz
                mdot_out[i] = mdot_out[i] - ro[k,j,i]*max(vx[k,j,i],0)*x[i]*dy*dz
                mdot_in[i]  = mdot_in[i]  - ro[k,j,i]*min(vx[k,j,i],0)*x[i]*dy*dz

    return mdot_net,mdot_out,mdot_in

def vertical_integrate(qq,depth,x,y,z,zpt=1):
    nx = qq.shape[2]
    ny = qq.shape[1]
    nz = qq.shape[0]

    sig = np.zeros([ny,nx])

    sf = surface_opt_av_azim(depth, z, zpt=zpt)
    ## if the tau=zpi dose not exist
    z0 = np.argmin(abs(z))
    for i in range(nx):
        if(abs(z[sf[0,i]]) < 0.5):
            sf[0,i] = np.argmin(abs(z+x[i]))
        if(abs(z[sf[1,i]]) < 0.5):
            sf[1,i] = np.argmin(abs(z-x[i]))
    
    for i in range(nx):
        for j in range(ny):
            for k in range(sf[0,i],sf[1,i]):
                dz = (z[k+1]-z[k-1])*0.5
                sig[j,i] = sig[j,i] + qq[k,j,i]*dz

    return sig

def cyl_surface_integrate(qq,depth,x,y,z,zpt=10):
    nx = qq.shape[2]
    ny = qq.shape[1]
    nz = qq.shape[0]

    sig = np.zeros([ny,nx])

    sf = surface_opt(depth, z, zpt=zpt)
   
    dy = abs(y[1]-y[0])
    for i in range(nx):
        for j in range(ny):
            for k in range(sf[0,j,i],sf[1,j,i]):
                dz = (z[k+1]-z[k-1])*0.5
                sig[j,i] = sig[j,i] + qq[k,j,i]*x[i]*dy*dz

    return sig

def volume_integrate(qq,depth,x,y,z,xi=0,xf=200,zpt=10):
    nx = qq.shape[2]
    ny = qq.shape[1]
    nz = qq.shape[0]

    z10= np.argmin(abs(z-zpt))
    z10m=np.argmin(abs(-z-zpt))
    z0 = np.argmin(abs(z))

    sf = surface_opt(depth, z, zpt=zpt)

    x1 = np.argmin(abs(x-xi))
    x2 = np.argmin(abs(x-xf))

    dy = y[1]-y[0]
    
    vol = 0.0
    for i in range(x1,x2):
        dx = (x[i+1]-x[i-1])
        for j in range(ny):
            for k in range(sf[0,j,i],sf[1,j,i]):
                dz = abs(z[k+1]-z[k-1])*0.5
                vol = vol + qq[k,j,i]*x[i]*dy*dx*dz

    return vol

#find optical depth = zpt
def surface_opt(tau, z, zpt=10):
    nx = tau.shape[2]
    ny = tau.shape[1]
    nz = tau.shape[0]

    surf = np.zeros([2,ny,nx])

    z0 = np.argmin(abs(z))

    for j in range(ny):
        for i in range(nx):
            tmp = tau[z0,j,i]-zpt
            surf[0,j,i] = np.argmin(abs(tau[:z0,j,i]-zpt)) 
            surf[0,j,i] = surf[0,j,i] + (z0+1 - surf[0,j,i])*min(tmp, 0)/tmp
            tmp = tau[z0+1,j,i]-zpt
            surf[1,j,i] = max(np.argmin(abs(tau[z0+1:,j,i]-zpt))+z0+1, z0*max(tmp, 0)/tmp)
    surf = np.asarray(surf,int)

    return surf


#find optical depth = zpt
def surface_opt_av_azim(tau, z, zpt=10):
    nx = tau.shape[2]
    nz = tau.shape[0]

    surf = np.zeros([2,nx])

    z0 = np.argmin(abs(z))
    mtau = np.average(tau,axis=1)

    for i in range(nx):
        tmp = mtau[z0,i]-zpt
        surf[0,i] = np.argmin(abs(mtau[:z0,i]-zpt)) 
        surf[0,i] = surf[0,i] + (z0+1 - surf[0,i])*min(tmp, 0)/tmp
        tmp = mtau[z0+1,i]-zpt
        surf[1,i] = max(np.argmin(abs(mtau[z0+1:,i]-zpt))+z0+1, z0*max(tmp, 0)/tmp)
    surf = np.asarray(surf,int)

    return surf
