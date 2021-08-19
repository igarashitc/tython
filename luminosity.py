import numpy as np

def lightcurve(fz,x,y,z,xf,zpt):
    nt = fz.shape[0]

    lmst = np.ndarray(nt*2).reshape(nt,2)

    for it in range(nt):
        lmst[it,:] = luminosity_cyl_z(fz[it,:,:,:],x,y,z,xf,zpt)

    Ledd = 1.25e45
    nmlf = 1e-10*27e30*9e24

    return lmst*nmlf/Ledd

#integrate out going z-direction flux in normalized unit
def luminosity_cyl_z(fz,x,y,z,xf,zpt):
    
    nx = fz.shape[2]
    ny = fz.shape[1]
    nz = fz.shape[0]

    lmst = np.ndarray(2)
    lmst[:] = 0

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
