
import numpy as np
import scipy as np
import matplotlib.pyplot as plt

class plot_cyl:

    def polar_2d(r, p, data, cmap="coolwarm", vmax=2, vmin=-2):
        dr = abs(r[1]-r[0])*0.5e0
        dp = abs(p[1]-p[0])*0.5e0

        tmp = p.copy()
        tmp[1:] = tmp[1:]+dp
        Rad,Phi = np.meshgrid(r-0.5*dr,tmp)
        X,Y = Rad*np.cos(Phi), Rad*np.sin(Phi)

        plt.pcolormesh(X,Y,data,cmap=cmap,vmax=vmax,vmin=vmin,shading="gourand")
        plt.gca().set_aspect("equal")

        return 

    #
    # plot poloidal plane
    #
    def plot_rz_spread(self, qql, qqr, x, z,\
                        x1=80., z1=-20., z2=60.,\
                        xts=np.linspace(0,80,5), zts=np.linspace(-20,60,5),\
                        cml="plasma", cmr="plasma",\
                        vmxl=0., vmnl=-5., vmxr=12., vmnr=6.,\
                        tsl=np.linspace(-5,0,6), tsr=np.linspace(6,12,7),\
                        ttl=r"$\log\ \rho/\rho_0$", ttr=r"$\log\ T\ \rm{[K]}$",\
                        time=1000):
    
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    
        fig, [ax1, ax2] = plt.subplots(1,2,figsize=[5,4])
        plt.subplots_adjust(wspace=0,left=0.15,right=0.85,bottom=0.1)
        plt.rcParams["font.size"]=12
    
        # left panel
        tmp = np.fliplr(qql).copy()
        xf  = -np.flip(x)
        axins1 = inset_axes(ax1, width="80%", height="5%", loc="lower left",\
                            bbox_transform=ax1.transAxes, bbox_to_anchor=(0,1.025, 1, 1))
        pc1 = ax1.pcolormesh(xf, z, tmp, cmap=cml, vmax=vmxl, vmin=vmnl)
        ax1.set_xlim(-x1,0)
        ax1.set_ylim(z1,z2)
        ax1.set_aspect("equal")
        ax1.set_xlabel(r"$r/r_{\rm s}$")
        ax1.set_ylabel(r"$z/r_{\rm s}$")
        ax1.set_xticks(-xts)
        ax1.set_yticks(zts)
        cb1 = fig.colorbar(pc1, cax=axins1, orientation="horizontal", ticks=tsl, label=ttl)
        axins1.xaxis.set_ticks_position("top")
        axins1.xaxis.set_label_position("top")
    
        # right panel
        axins2 = inset_axes(ax2, width="80%", height="5%", loc="lower left",\
                            bbox_transform=ax2.transAxes, bbox_to_anchor=(0.125, 1.025, 1, 1))
        pc2 = ax2.pcolormesh(x, z, qqr, cmap=cmr, vmax=vmxr, vmin=vmnr)
        ax2.set_xlim(0,x1)
        ax2.set_ylim(z1,z2)
        ax2.set_aspect("equal")
        ax2.tick_params(labelleft=False, labelright=True)
        ax2.set_xlabel(r"$r/r_{\rm s}$")
        #ax2.set_ylabel(r"$z/r_{\rm s}$")
        ax2.set_xticks(xts)
        ax2.set_yticks(zts)
        cb2 = fig.colorbar(pc2, cax=axins2, orientation="horizontal",ticks=tsr, label=ttr)
        axins2.xaxis.set_ticks_position("top")
        axins2.xaxis.set_label_position("top")
    
        #time
        #fig.suptitle(str(time)+r"$t_0$")
        fig.suptitle(time)
    
        return fig
    
    def movie_rz(self,data,x,y,z,cmap,vmx,vmn,x1,z1,z2,dir):
        nt = data.shape[0]
    
        for l in range(0,nt-1):
            plt.pcolormesh(x,z,data[l,:,0,:],cmap=cmap,vmax=vmx,vmin=vmn)
            plt.xlim(0,x1)
            plt.ylim(z1,z2)
            plt.gca().set_aspect("equal",adjustable='box')
            plt.colorbar()
            plt.savefig(dir+str(l).zfill(3)+'.png',bbox_inches='tight',dpi=150)
            plt.close()
    
        return
    
    # plot psd at each radius.
    def psd_scipy_2d(self,data,time,dt):
        import scipy.fftpack as sf
    
        nx   = data.shape[1]
    
        N    = len(time)
        posN = len(time[0:int(N/2)])
    
        psd  = np.zeros([posN, nx])
    
        freq = sf.fftfreq(N, dt)[0:int(N/2)]
        df   = freq[1]-freq[0]
    
        for i in range(nx):
            fft      = sf.fft(data[:,i], N)[0:int(N/2)]
            psd[:,i] = np.abs(fft)/np.sqrt(df)/(N/np.sqrt(2))
    
        return psd, freq
