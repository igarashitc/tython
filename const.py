# physical constants in cgs units

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
#unit density in this simulation
de = 1e-10
#black hole mass per solar mass
bh = 1e7
#proton mass 
mp = 1.6726218e-24
#electron mass 
me = 0
#specific head ratio
gm = 5./3.
#gravitational constant
gc = 6.67428e-8
#shwarzshild radius
rs = 3e5*bh
#solar mass
Msun = 1.98841586057e33

#Eddington luminosity
Ledd = 4*np.pi*cc*gc*bh*Msun/0.4
#Eddington accretionrate 
mded = Ledd/cc**2

#normalized temperature in this simulation
te0 = mmw*mp*cc**2/kb
    
