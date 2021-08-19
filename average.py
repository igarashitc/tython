import numpy as np

def azimuthal(qq,x,y,z):
    nx = qq.shape[2]
    ny = qq.shape[1]
    nz = qq.shape[0]

    mqq = np.ndarray(nz*nx).reshape(nz,nx)

    for k in range(nz):
        for i in range(nx):
            mqq[k,i] = np.average(qq[k,:,i]) 

    return mqq
