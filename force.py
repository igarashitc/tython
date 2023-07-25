import numpy as np

def calc(ro,pr,vx,vy,vz,bx,by,bz,fx,fy,fz,er,x,y,z,density=1,length=1):
    from tython import radiation as rd

    nx = ro.shape[2]
    ny = ro.shape[1]
    nz = ro.shape[0]
    tmp = np.ones([nz,ny,nx])

    #innertial term (drT_rr/rdr + dT_phir/rdphi + dT_zr/dz - ro*v_phi^2/r) 
    mom = -adv_cyl_tensor_r(vx,vy,vz,x,y,z)
    #magnetic tension 
    mag = div_cyl_tensor_r(bx,by,bz,x,y,z)
    mag = mag
    #gas pressure gradient (dp/dr)
    pgf = -grad_r(pr,x)
    pb  = (bx**2+by**2+bz**2)*0.5
    #magnetic pressure gradient (dp/dr)
    pbf = -grad_r(pb,x)
    #gravitational force
    grf = -grav_r(tmp,x,z)
    #radiation forcd
    rad = -rd.force_r(ro,pr,vx,vy,vz,er,fx,fy,fz,density=density,length=length) 
    #centrifugal force (v_phi^2/r)
    ctf = cent(vy,x)
    #magnetic centrifugal force
    mcf = -cent(by,x)

    for k in range(nz-1):
        for j in range(ny-1):
            for i in range(1,nx-1):
                roi        = 1/ro[k,j,i]
                mag[k,j,i] = mag[k,j,i]*roi
                mcf[k,j,i] = mcf[k,j,i]*roi
                pgf[k,j,i] = pgf[k,j,i]*roi
                pbf[k,j,i] = pbf[k,j,i]*roi

    tot = mom+mag+pgf+pbf+grf+rad+ctf+mcf

    return mom,mag,pgf,pbf,grf,rad,ctf,mcf,tot


def cent(qq,x):
    nx = qq.shape[2]
    ny = qq.shape[1]
    nz = qq.shape[0]

    cent = np.ndarray(nz*ny*nx).reshape(nz,ny,nx)

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cent[k,j,i] = qq[k,j,i]*qq[k,j,i]/x[i]

    return cent

def grav_r(qq,x,z):
    nx = qq.shape[2]
    ny = qq.shape[1]
    nz = qq.shape[0]

    grav = np.ndarray(nz*ny*nx).reshape(nz,ny,nx)
    
    #pseudo-Newtonian potential
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                rr = np.sqrt(x[i]**2+z[k]**2)
                grav[k,j,i] = (qq[k,j,i]*x[i] / (2*(rr-1)**2*rr))

    return grav

# dqq/dr
def grad_r(qq,x):
    nx = qq.shape[2]
    ny = qq.shape[1]
    nz = qq.shape[0]

    grad = np.ndarray(nz*ny*nx).reshape(nz,ny,nx)

    for k in range(nz-1):
        for j in range(ny-1):
            for i in range(1,nx-1):
                dx =1./(x[i+1] - x[i-1]) 
                grad[k,j,i] = (qq[k,j,i+1] - qq[k,j,i-1])*dx

    return grad

# d/rdr (rT_rr) + d/rdphi (T_phir) + d/dz (T_zr)
def div_cyl_tensor_r(qx,qy,qz,x,y,z):

    nx = qx.shape[2]
    ny = qx.shape[1]
    nz = qx.shape[0]

    divq = np.ndarray(nz*ny*nx).reshape(nz,ny,nx)

    pqx = periodic(qx)
    pqy = periodic(qy)

    dy = 1/np.abs(y[2]-y[0])

    for k in range(1,nz-1):
        dz = 1./(z[k+1]-z[k-1])
        for j in range(0,ny-1):
            for i in range(1,nx-1):
                dx = 1/(x[i+1]-x[i-1])
                xi = 1/x[i]
                dqdx = (x[i+1]*qx[k,j,i+1]**2 - x[i-1]*qx[k,j,i-1]**2)*dx*xi
                dqdy = (pqx[k,j+2,i]*pqy[k,j+2,i] - pqx[k,j,i]*pqy[k,j,i])*dy*xi
                dqdz = (qx[k+1,j,i]*qz[k+1,j,i] - qx[k-1,j,i]*qz[k-1,j,i])*dz
                divq[k,j,i] = dqdx+dqdy+dqdz

    return divq

# v_r*(dv_r/dr) + v_phi/r*(dv_r/dphi) + v_z*dv_r/dz (q (\cdot\nabla) q)
def adv_cyl_tensor_r(qx,qy,qz,x,y,z):

    nx = qx.shape[2]
    ny = qx.shape[1]
    nz = qx.shape[0]

    divq = np.zeros([nz,ny,nx])

    pqx = periodic(qx)
    pqy = periodic(qy)

    dy = 1/np.abs(y[2]-y[0])

    for k in range(1,nz-1):
        dz = 1./(z[k+1]-z[k-1])
        for j in range(0,ny-1):
            for i in range(1,nx-1):
                dx = 1/(x[i+1]-x[i-1])
                xi = 1/x[i]
                dqdx = qx[k,j,i]*(qx[k,j,i+1] - qx[k,j,i-1])*dx
                dqdy = (pqx[k,j+2,i]*pqy[k,j+2,i] - pqx[k,j,i]*pqy[k,j,i])*dy*xi
                dqdy = qy[k,j,i]*xi*(pqx[k,j+2,i] - pqx[k,j,i])*dy*xi
                dqdz = qz[k,j,i]*(qx[k+1,j,i] - qx[k-1,j,i])*dz
                divq[k,j,i] = dqdx+dqdy+dqdz

    return divq

def rotation_cyl_r(qx,qy,qz,x,y,z):

    nx = qx.shape[2]
    ny = qx.shape[1]
    nz = qx.shape[0]

    rotq = np.zeros(nz,ny,nx)

    pqz = periodic(qz)

    dy = 1/np.abs(y[2]-y[0])

    for k in range(1,nz-1):
        dz = 1./(z[k+1]-z[k-1])
        for j in range(0,ny-1):
            for i in range(1,nx-1):
                xi = 1/x[i]
                rotq[k,j,i] = (pqz[k,j+2,i] - pqz[k,j,i])*dy*xi
                rotq[k,j,i] = rotq[k,j,i] - (qy[k+1,j,i] - qy[k,j,i])*dz

    return rotq

def periodic(qq):
    nx = qq.shape[2]
    ny = qq.shape[1]
    nz = qq.shape[0]

    pqq = np.ndarray(nz*(ny+2)*nx).reshape(nz,ny+2,nx)

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                pqq[k,j+1,i] = qq[k,j,i]

    pqq[:,0   ,:] = qq[:,ny-1,:]
    pqq[:,ny+1,:] = qq[:,0   ,:]

    return pqq

