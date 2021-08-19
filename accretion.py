import numpy as np

def calc(ro,vx,x,y,z,density,length):
    from python import radiation as rd

    nt = ro.shape[0]
    nx = ro.shape[3]

    mdot = np.ndarray(nt*nx).reshape(nt,nx)
    sig  = np.ndarray(nt*nx).reshape(nt,nx)
    
    for it in range(nt):
        tde = rd.thomson_depth(ro[it,:,:,:],z,density,length)
        mdot[it,:] = mdot_net(ro[it,:,:,:],vx[it,:,:,:],tde,x,y,z)
        sig[it,:] = surface_density(ro[it,:,:,:],tde,x,y,z)
        print(it)

    return mdot,sig

def mdot_net(ro,vx,depth,x,y,z):
    nx = ro.shape[2]
    ny = ro.shape[1]
    nz = ro.shape[0]

    mdot = np.ndarray(nx)
    mdot[:] = 0.

    dy = abs(y[1]-y[0])

    zpt= 10.
    z10= np.argmin(abs(z-10))
    z10m=np.argmin(abs(-z-10))
    z0 = np.argmin(abs(z))

    sf = np.ndarray(2*ny*nx).reshape(2,ny,nx)
    for j in range(ny):
        for i in range(nx):
            sf[0,j,i] = min(z10m, np.argmin(abs(depth[:z0+1,j,i]-10)))
            sf[1,j,i] = max(z10, np.argmin(abs(depth[z0+1:,j,i]-10)) + z0+1)
    sf = np.asarray(sf,int)

    for i in range(nx-1):
        for j in range(ny):
            for k in range(sf[0,j,i],sf[1,j,i]):
                dz = (z[k+1]-z[k-1])*0.5
                mdot[i] = mdot[i] -ro[k,j,i]*vx[k,j,i]*x[i]*dy*dz

    return mdot

def surface_density(ro,depth,x,y,z):
    nx = ro.shape[2]
    ny = ro.shape[1]
    nz = ro.shape[0]

    sig = np.ndarray(nx)
    sig[:] = 0.

    zpt= 40.
    zu = np.argmin(abs(z-zpt))
    zd = np.argmin(abs(-z-zpt))

    zpt= 10.
    z10= np.argmin(abs(z-10))
    z10m=np.argmin(abs(-z-10))
    z0 = np.argmin(abs(z))

    sf = np.ndarray(2*ny*nx).reshape(2,ny,nx)
    for j in range(ny):
        for i in range(nx):
            sf[0,j,i] = min(z10m, np.argmin(abs(depth[:z0+1,j,i]-10)))
            sf[1,j,i] = max(z10, np.argmin(abs(depth[z0+1:,j,i]-10)) + z0+1)
    sf = np.asarray(sf,int)

    for k in range(zd,zu):
        dz = (z[k+1]-z[k-1])*0.5
        for j in range(ny):
            for i in range(nx):
                sig[i] = sig[i] + ro[k,j,i]*dz

    return sig
