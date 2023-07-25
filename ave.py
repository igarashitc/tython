
t1 = 10
t2 = 25
t3 = 79
t4 = 85
t5 = 125
t6 = 139

ro = np.load("npy/ro.npy")
tro1 = average(ro[t1:t2,:,:,:],axis=0)
tro2 = average(ro[t3:t4,:,:,:],axis=0)
tro3 = average(ro[t5:t6,:,:,:],axis=0)

vx = np.load("npy/vx.npy")
tvx1 = average(vx[t1:t2,:,:,:],axis=0)
tvx2 = average(vx[t3:t4,:,:,:],axis=0)
tvx3 = average(vx[t5:t6,:,:,:],axis=0)

#del ro,vx
del vx
vy = np.load("npy/vy.npy")
tvy1 = average(vy[t1:t2,:,:,:],axis=0)
tvy2 = average(vy[t3:t4,:,:,:],axis=0)
tvy3 = average(vy[t5:t6,:,:,:],axis=0)

vz = np.load("npy/vz.npy")
tvz1 = average(vz[t1:t2,:,:,:],axis=0)
tvz2 = average(vz[t3:t4,:,:,:],axis=0)
tvz3 = average(vz[t5:t6,:,:,:],axis=0)

#del vy,vz
del vz
bx = np.load("npy/bx.npy")
tbx1 = average(bx[t1:t2,:,:,:],axis=0)
tbx2 = average(bx[t3:t4,:,:,:],axis=0)
tbx3 = average(bx[t5:t6,:,:,:],axis=0)

by = np.load("npy/by.npy")
tby1 = average(by[t1:t2,:,:,:],axis=0)
tby2 = average(by[t3:t4,:,:,:],axis=0)
tby3 = average(by[t5:t6,:,:,:],axis=0)
#del bx,by

bz = np.load("npy/bz.npy")
tbz1 = average(bz[t1:t2,:,:,:],axis=0)
tbz2 = average(bz[t3:t4,:,:,:],axis=0)
tbz3 = average(bz[t5:t6,:,:,:],axis=0)

er = np.load("npy/er.npy")
ter1 = average(er[t1:t2,:,:,:],axis=0)
ter2 = average(er[t3:t4,:,:,:],axis=0)
ter3 = average(er[t5:t6,:,:,:],axis=0)
#del bz,er

pr = np.load("npy/pr.npy")
tpr1 = average(pr[t1:t2,:,:,:],axis=0)
tpr2 = average(pr[t3:t4,:,:,:],axis=0)
tpr3 = average(pr[t5:t6,:,:,:],axis=0)

fx = np.load("npy/fx.npy")
tfx1 = average(fx[t1:t2,:,:,:],axis=0)
tfx2 = average(fx[t3:t4,:,:,:],axis=0)
tfx3 = average(fx[t5:t6,:,:,:],axis=0)
#del pr,fx
del fx

fy = np.load("npy/fy.npy")
tfy1 = average(fy[t1:t2,:,:,:],axis=0)
tfy2 = average(fy[t3:t4,:,:,:],axis=0)
tfy3 = average(fy[t5:t6,:,:,:],axis=0)

fz = np.load("npy/fz.npy")
tfz1 = average(fz[t1:t2,:,:,:],axis=0)
tfz2 = average(fz[t3:t4,:,:,:],axis=0)
tfz3 = average(fz[t5:t6,:,:,:],axis=0)
del fy,fz
