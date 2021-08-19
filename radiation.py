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


def force_r(ro,pr,vx,vy,vz,er,fx,fy,fz,density,length):
    nx = ro.shape[2]
    ny = ro.shape[1]
    nz = ro.shape[0]

    #normalized radiation constant
    ar0 = ar*mmw**4*mp**4*cc**6/(density*kb**4)
    #normalized electron scattering coeff
    kes = 0.4*density*length

    kap = absorption(ro,pr,density,length)

    sr = np.ndarray(nz*ny*nx).reshape(nz,ny,nx)

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

                tmp1 = ro[k,j,i]*kap[k,j,i]*vx[k,j,i]*(ar0*(pr[k,j,i]/ro[k,j,i])**4 - er[k,j,i])
                tmp2 = ro[k,j,i]*(kap[k,j,i]+kes)
                tmp2 = tmp2*(fx[k,j,i] - (vx[k,j,i]*er[k,j,i] + vx[k,j,i]*prr+vy[k,j,i]*prp+vz[k,j,i]*prz))

                sr[k,j,i] = tmp1 - tmp2

    return sr

def Eddington_tensor_m1(er,fx,fy,fz)
    
    Dr = np.ndarray(3*3).reshape(3,3)

    f2 = fx**2+fy**2+fz**2
    ff = f2/er[k,j,i]**2
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

def absorption(ro,pr,density,length):
    kp0 = 0.64e23*density*length
    
    absorption = kp0*density*(kb/(mmw*mp*cc*cc))**3.5
    absorption = absorption*ro*(pr/ro)**(-3.5)
    
    return absorption

def thomson_depth(ro,z,density,length):
    nx = ro.shape[2]
    ny = ro.shape[1]
    nz = ro.shape[0]

    kes = 0.4*density*length

    depth = np.ndarray(nz*ny*nx).reshape(nz,ny,nx)
    depth[:,:,:] = 0.

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
