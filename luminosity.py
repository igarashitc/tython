import numpy as np

cc = 3e10

def lightcurve(fz,x,y,z,xf=200,zpt=100,Mbh=10,density=1,length=1):
    nt = fz.shape[0]

    lmst = np.zeros([nt,2])

    for it in range(nt):
        lmst[it,:] = luminosity_cyl_z(fz[it,:,:,:],x,y,z,xf=xf,zpt=zpt)

    Ledd = 1.25e38*Mbh
    nmlf = density*cc**3*length**2

    return lmst*nmlf/Ledd

#integrate out going z-direction flux in normalized unit
def luminosity_cyl_z(fz,x,y,z,xf=200,zpt=100):
    
    nx = fz.shape[2]
    ny = fz.shape[1]
    nz = fz.shape[0]

    print(zpt)

    lmst = np.zeros(2)

    dy = np.abs(y[1]-y[0])

    x1 = np.argmin(np.abs(x-xf))
    zu = np.argmin(np.abs(z-zpt))
    zd = np.argmin(np.abs(-z-zpt))

    for j in range(ny):
        for i in range(0,x1):
            dx = np.abs(x[i+1]-x[i-1])
            lmst[0] = lmst[0] + max(0,fz[zu,j,i])*0.5*dx*x[i]*dy
            lmst[1] = lmst[1] + np.abs(min(0,fz[zd,j,i]))*0.5*dx*x[i]*dy


    return lmst

#integrate out going z-direction flux in normalized unit
def luminosity_cyl_tau(fz,tde,x,y,z,xi=0,xf=200,zpt=1,Mbh=1e7,density=1,length=1):

    from tython import accretion as ac
    
    nx = fz.shape[2]
    ny = fz.shape[1]
    nz = fz.shape[0]

    lmst = np.zeros(2)

    dy = np.abs(y[1]-y[0])

    x1 = np.argmin(abs(x-xi))
    x2 = np.argmin(abs(x-xf))
    z0 = np.argmin(abs(z))

    sf = ac.surface_opt(tde,z,zpt=zpt)

    for j in range(ny):
        for i in range(x1,x2):
            dx = np.abs(x[i+1]-x[i-1])
            lmst[0] = lmst[0] + max(0,fz[sf[0,j,i],j,i])*0.5*dx*x[i]*dy
            lmst[1] = lmst[1] + abs(min(0,fz[sf[1,j,i],j,i]))*0.5*dx*x[i]*dy

    Ledd = 1.25e38*Mbh
    nmlf = density*cc**3*length**2

    return lmst*nmlf/Ledd
