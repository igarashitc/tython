import numpy as np

# volume average
def cyl_volume(qq,x,y,z,x1=0,x2=40,z1=-10,z2=10):

    ny = qq.shape[1]

    dy = np.abs(y[1]-y[0])

    mqq = 0.
    vol = 0.
    for k in range(np.argmin(abs(x-x1)), np.argmin(abs(x-x2))):
        dz = np.abs(z[k+1]-z[k-1])*0.5
        for j in range(ny):
            for i in range(np.argmin(abs(z-z1)), np.argmin(abs(z-z2))):
                dx  = np.abs(x[i+1]-x[i-1])*0.5
                tmp = x[i]*dy*dx*dz
                vol = vol + tmp
                mqq = mqq + qq[k,j,i]*tmp 

    return mqq/vol

# volume average within an optical depth
def cyl_volume_opt(qq,depth,x,y,z,x1=0,x2=40,zpt=1):
    from tython import accretion 

    nx = qq.shape[2]
    ny = qq.shape[1]

    dy = np.abs(y[1]-y[0])

    xs = np.argmin(abs(x-x1))
    xf = np.argmin(abs(x-x2))

    sf = accretion.surface_opt(depth,z,zpt=zpt)

    mqq = 0.
    vol = 0.
    for j in range(ny):
        for i in range(xs, xf):
            dx  = np.abs(x[i+1]-x[i-1])*0.5
            for k in range(sf[0,j,i], sf[1,j,i]):
                dz = np.abs(z[k+1]-z[k-1])*0.5

                tmp = x[i]*dy*dx*dz
                vol = vol + tmp
                mqq = mqq + qq[k,j,i]*tmp 

    print(vol)
    return mqq/vol

# volume average within an optical depth
def cyl_surface_opt(qq,depth,x,y,z,zpt=1):
    from tython import accretion 

    nx = qq.shape[2]
    ny = qq.shape[1]

    dy = np.abs(y[1]-y[0])

    sf = accretion.surface_opt(depth,z,zpt=zpt)

    mqq = np.zeros(nx)
    vol = np.zeros(nx)
    for j in range(ny):
        for i in range(0, nx-1):
            for k in range(sf[0,j,i], sf[1,j,i]):
                dz = np.abs(z[k+1]-z[k-1])*0.5

                tmp = dy*dz
                vol[i] = vol[i] + tmp
                mqq[i] = mqq[i] + qq[k,j,i]*tmp 

    return mqq/vol,vol
