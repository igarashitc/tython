import numpy
from scipy.interpolate import RegularGridInterpolator as rgi

def cyl2cart_3d(dat_cyl,r_cyl,phi_cyl,z_cyl):

    nr = dat_cyl.shape[2]
    np = dat_cyl.shape[1]
    nz = dat_cyl.shape[0]

    tmp           = numpy.zeros([nz,np+1,nr])
    tmp[:,0:np,:] = dat_cyl
    tmp[:,np  ,:] = dat_cyl[:,0,:]

    nx = nr*2
    ny = nr*2
    x  = numpy.ndarray(nx)
    y  = numpy.ndarray(ny)
    z  = numpy.ndarray(nz)
    dat= numpy.zeros([nz,ny,nx])

    rmax = numpy.max(r_cyl) 
    rmin = numpy.min(r_cyl)
    zmax = numpy.max(z_cyl)
    zmin = numpy.min(z_cyl)
    x[nr:nx  ]  = 0.0 + (numpy.arange(nr)+0.5)*(rmax-rmin)/(nr-1.)
    x[0:nr   ]  = -numpy.flipud(x[nr:nx])
    y = x
    z = zmin+numpy.arange(nz)*(zmax-zmin)/(nz-1.)

    r_tmp           = numpy.ndarray(nr+1)
    r_tmp[0:nr  ]   = r_cyl
    r_tmp[nr    ]   = 2*r_cyl[nr-1]-r_cyl[nr-2]
    phi_tmp         = numpy.ndarray(np+1)
    phi_tmp[0:np]   = phi_cyl
    phi_tmp[np  ]   = 2*phi_cyl[np-1]-phi_cyl[np-2]
    z_tmp           = numpy.ndarray(nz+1)
    z_tmp[0:nz  ]   = z_cyl
    z_tmp[nz    ]   = 2*z_cyl[nz-1]-z_cyl[nz-2]

    nzs = 0
    nze = nz

    fn = rgi((z_cyl,phi_tmp,r_cyl), tmp)
    for k in range(nzs,nze):
        print(k)
        for j in range(0,ny):
            sgn = numpy.sign(y[j])
            for i in range(0,nx):
                r   = numpy.sqrt(x[i]**2 + y[j]**2)
                phi = numpy.arccos(x[i]*sgn/r) - numpy.pi*numpy.min([0,sgn])
                r   = min([r,rmax])

                dat[k,j,i] = fn([z[k], phi, r])


    return dat,x,y,z

def cyl2cart_3d_mid(dat_cyl,r_cyl,phi_cyl):

    nr = dat_cyl.shape[1]
    np = dat_cyl.shape[0]

    tmp           = numpy.zeros([np+1,nr])
    tmp[0:np,:] = dat_cyl
    tmp[np  ,:] = dat_cyl[0,:]

    nx = nr*2
    ny = nr*2
    x  = numpy.ndarray(nx)
    y  = numpy.ndarray(ny)
    dat= numpy.zeros([ny,nx])

    rmax = numpy.max(r_cyl) 
    rmin = numpy.min(r_cyl)
    x[nr:nx  ]  = 0.0 + (numpy.arange(nr)+0.5)*(rmax-rmin)/(nr-1.)
    x[0:nr   ]  = -numpy.flipud(x[nr:nx])
    y = x

    r_tmp           = numpy.ndarray(nr+1)
    r_tmp[0:nr  ]   = r_cyl
    r_tmp[nr    ]   = 2*r_cyl[nr-1]-r_cyl[nr-2]
    phi_tmp         = numpy.ndarray(np+1)
    phi_tmp[0:np]   = phi_cyl
    phi_tmp[np  ]   = 2*phi_cyl[np-1]-phi_cyl[np-2]

    fn = rgi((phi_tmp,r_cyl), tmp)
    for j in range(0,ny):
        sgn = numpy.sign(y[j])
        for i in range(0,nx):
            r   = numpy.sqrt(x[i]**2 + y[j]**2)
            phi = numpy.arccos(x[i]*sgn/r) - numpy.pi*numpy.min([0,sgn])
            r   = min([r,rmax])

            dat[j,i] = fn([phi, r])

    return dat,x,y

def cyl2cart_2d(dat_cyl,r_cyl,z_cyl):

    nr = dat_cyl.shape[1]
    nz = dat_cyl.shape[0]

    tmp      = numpy.zeros([nz,nr])
    tmp[:,:] = dat_cyl
    tmp[:,:] = dat_cyl[:,:]

    nx = nr*2
    x  = numpy.ndarray(nx)
    z  = numpy.ndarray(nz)
    dat= numpy.zeros([nz,nx])

    rmax = numpy.max(r_cyl) 
    rmin = numpy.min(r_cyl)
    zmax = numpy.max(z_cyl)
    zmin = numpy.min(z_cyl)
    x[nr:nx  ]  = 0.0 + (numpy.arange(nr)+0.5)*(rmax-rmin)/(nr-1.)
    x[0:nr   ]  = -numpy.flipud(x[nr:nx])
    z = zmin+numpy.arange(nz)*(zmax-zmin)/(nz-1.)

    r_tmp           = numpy.ndarray(nr+1)
    r_tmp[0:nr  ]   = r_cyl
    r_tmp[nr    ]   = 2*r_cyl[nr-1]-r_cyl[nr-2]
    z_tmp           = numpy.ndarray(nz+1)
    z_tmp[0:nz  ]   = z_cyl
    z_tmp[nz    ]   = 2*z_cyl[nz-1]-z_cyl[nz-2]

    nzs = 0
    nze = nz

    fn = rgi((z_cyl,r_cyl), tmp)
    for k in range(nzs,nze):
        for i in range(0,nx):
            r   = numpy.sqrt(x[i]**2)
            r   = min([r,rmax])

            dat[k,i] = fn([z[k], r])


    return dat,x,z
