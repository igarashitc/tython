import numpy as np
import matplotlib.pyplot as plt

def movie(data,x,y,z,cmap,vmx,vmn,x1,z1,z2,dir):
    nt = data.shape[0]

    for l in range(0,nt):
        plt.pcolormesh(x,z,data[l,:,0,:],cmap=cmap,vmax=vmx,vmin=vmn)
        plt.xlim(0,x1)
        plt.ylim(z1,z2)
        plt.gca().set_aspect("equal",adjustable='box')
        plt.colorbar()
        plt.savefig(dir+str(l).zfill(3)+'.png',bbox_inches='tight',dpi=150)
        plt.close()


