import numpy as np

#stefan-boltzman constant
sb = 5.6703730e-5
#speed of light 
cc = 3e10
#boltzmann constatn
kb = 1.3806488e-16
#radiation constatn
ar = 4*sb/cc
#mean molecular weight
mmw = 0.5
#proton mass 
mp = 1.6726218e-24
#electron mass 
me = 9.10938188e-28

def force_r(ro,pr,vx,vy,vz,er,fx,fy,fz,density=1,length=1):
    nx = ro.shape[2]
    ny = ro.shape[1]
    nz = ro.shape[0]

    #normalized radiation constant
    ar0 = ar*mmw**4*mp**4*cc**6/(density*kb**4)
    #normalized electron scattering coeff
    kes = 0.4*density*length

    kap = ff_absorption(ro,pr,density=density,length=length)

    sr = np.zeros([nz,ny,nx])

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                f2 = fx[k,j,i]**2+fy[k,j,i]**2+fz[k,j,i]**2
                ff = f2/er[k,j,i]**2
                chi= (3 + 4*ff)/(5 + 2*np.sqrt(4-3*ff))
                f2 = np.sqrt(f2)
                lx = fx[k,j,i]/f2
                ly = fy[k,j,i]/f2
                lz = fz[k,j,i]/f2

                prr = ((1-chi)*0.5 + (3*chi-1)*0.5*lx*lx)*er[k,j,i]
                prp = ((3*chi-1)*0.5*lx*ly)*er[k,j,i]
                prz = ((3*chi-1)*0.5*lx*lz)*er[k,j,i]

                tmp1 = kap[k,j,i]*vx[k,j,i]*(ar0*(pr[k,j,i]/ro[k,j,i])**4 - er[k,j,i])
                tmp2 = kap[k,j,i]+kes
                tmp2 = tmp2*(fx[k,j,i] - (vx[k,j,i]*er[k,j,i] + vx[k,j,i]*prr+vy[k,j,i]*prp+vz[k,j,i]*prz))

                sr[k,j,i] = tmp1 - tmp2

    return sr

def Eddington_tensor_m1(er,fx,fy,fz):
    
    Dr = np.ndarray(3*3).reshape(3,3)

    f2 = fx**2+fy**2+fz**2
    ff = f2/er**2
    chi= (3 + 4*ff)/(5 + 2*np.sqrt(4-3*ff))
    f2 = np.sqrt(f2)
    lx = fx/f2
    ly = fy/f2
    lz = fz/f2

    Ddiag = 0.5*(1. - chi)
    Doff = 0.5*(3.*chi - 1.)
    Dr[0,0] = Doff*lx*lx + Ddiag
    Dr[1,0] = Doff*lx*ly
    Dr[2,0] = Doff*lx*lz
    Dr[0,1] = Dr[1,0]
    Dr[1,1] = Doff*ly*ly + Ddiag
    Dr[2,1] = Doff*lz*ly
    Dr[0,2] = Dr[2,0]
    Dr[1,2] = Dr[2,1]
    Dr[2,2] = Doff*lz*lz + Ddiag

    return Dr

def Eddington_tensor_m1_ar(er,fx,fy,fz):

    nx = er.shape[2]
    ny = er.shape[1]
    nz = er.shape[0]

    Dr = np.zeros([3,3,nz,ny,nx])

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                f2 = fx[k,j,i]**2+fy[k,j,i]**2+fz[k,j,i]**2
                ff = f2/er[k,j,i]**2
                chi= (3 + 4*ff)/(5 + 2*np.sqrt(4-3*ff))
                f2 = np.sqrt(f2)
                lx = fx[k,j,i]/f2
                ly = fy[k,j,i]/f2
                lz = fz[k,j,i]/f2

                Ddiag = 0.5*(1. - chi)
                Doff = 0.5*(3.*chi - 1.)
                Dr[0,0,k,j,i] = Doff*lx*lx + Ddiag
                Dr[1,0,k,j,i] = Doff*lx*ly
                Dr[2,0,k,j,i] = Doff*lx*lz
                Dr[0,1,k,j,i] = Dr[1,0,k,j,i]
                Dr[1,1,k,j,i] = Doff*ly*ly + Ddiag
                Dr[2,1,k,j,i] = Doff*lz*ly
                Dr[0,2,k,j,i] = Dr[2,0,k,j,i]
                Dr[1,2,k,j,i] = Dr[2,1,k,j,i]
                Dr[2,2,k,j,i] = Doff*lz*lz + Ddiag

    return Dr

def ff_absorption(ro,pr,density=1,length=1):
    kp0 = 0.64e23*density*length
    
    absorption = kp0*density*(kb/(mmw*mp*cc*cc))**3.5
    absorption = absorption*ro*(pr/ro)**(-3.5)
    
    return absorption

def absorption_depth(ro,pr,z,density=1,length=1):
    nx = ro.shape[2]
    ny = ro.shape[1]
    nz = ro.shape[0]

    kp0 = ff_absorption(ro,pr,density=density,length=length)

    #depth = np.zeros([nz,ny,nx])
    depth = np.zeros([nz,ny,nx])

    zpt = 100
    zu  = np.argmin(abs(z-zpt))
    zd  = np.argmin(abs(-z-zpt))
    z0  = np.argmin(abs(z))

    for k in range(zd,z0):
        dz = np.abs(z[k+1]-z[k-1])*0.5
        for j in range(ny):
            for i in range(nx):
                depth[k,j,i] = depth[k-1,j,i] + ro[k,j,i]*kp0[k,j,i]*dz

    for k in reversed(range(z0,zu)):
        dz = abs(z[k+1]-z[k-1])*0.5
        for j in range(ny):
            for i in range(nx):
                depth[k,j,i] = depth[k+1,j,i] + ro[k,j,i]*kp0[k,j,i]*dz

    return depth

def thomson_depth(ro,z,density=1,length=1):
    nx = ro.shape[2]
    ny = ro.shape[1]
    nz = ro.shape[0]

    kes = 0.4*density*length

    depth = np.zeros([nz,ny,nx])

    zpt = 100
    zu  = np.argmin(abs(z-zpt))
    zd  = np.argmin(abs(-z-zpt))
    z0  = np.argmin(abs(z))

    for k in range(zd,z0):
        dz = np.abs(z[k+1]-z[k-1])*0.5
        for j in range(ny):
            for i in range(nx):
                depth[k,j,i] = depth[k-1,j,i] + ro[k,j,i]*kes*dz

    for k in reversed(range(z0,zu)):
        dz = abs(z[k+1]-z[k-1])*0.5
        for j in range(ny):
            for i in range(nx):
                depth[k,j,i] = depth[k+1,j,i] + ro[k,j,i]*kes*dz

    return depth

def effective_depth(ro, pr, z, density=1, length=1):
    nx = ro.shape[2]
    ny = ro.shape[1]
    nz = ro.shape[0]

    depth = np.zeros([nz,ny,nx])    

    zpt = 100
    zu  = np.argmin(abs(z-zpt))
    zd  = np.argmin(abs(-z-zpt))
    z0  = np.argmin(abs(z))

    kab = ff_absorption(ro,pr,density=density,length=length)
    kes = 0.4*density*length

    for k in range(zd,z0):
        dz = np.abs(z[k+1]-z[k-1])*0.5
        for j in range(ny):
            for i in range(nx):
                depth[k,j,i] = depth[k-1,j,i] + ro[k,j,i]*dz*np.sqrt(kab[k,j,i]*(kab[k,j,i]+kes) )

    for k in reversed(range(z0,zu)):
        dz = abs(z[k+1]-z[k-1])*0.5
        for j in range(ny):
            for i in range(nx):
                depth[k,j,i] = depth[k+1,j,i] + ro[k,j,i]*dz*np.sqrt(kab[k,j,i]*(kab[k,j,i]+kes))

    return depth

def ic_cooling(ro,pr,er,density=1,length=1):

    #normalized radiation constant
    ar0 = ar*mmw**4*mp**4*cc**6/(density*kb**4)
    #normalized electron scattering coeff
    kes = 0.4*density*length
    
    cof= mmw*mp/me

    tr  = (er/ar0)**0.25
    tg  = pr/ro
    ic  = ro*kes*cof*er*(tg-tr)

    return ic

def timescale_ic(ro,pr,er,density=1,length=1):
   
    #normalized radiation constant
    ar0 = ar*mmw**4*mp**4*cc**6/(density*kb**4)
    #normalized electron scattering coeff
    kes = 0.4*density*length
    
    gm = 5./3.
    cof= 4*mmw*mp/me

    ee  = pr/(gm-1)
    tr  = (er/ar0)**0.25
    tg  = pr/ro
    ic  = ro*kes*cof*er*(tg-tr)
    tic = ee/ic

    return tic

