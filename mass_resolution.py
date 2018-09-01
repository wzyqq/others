import numpy as np

omega_m  = 0.268
boxsize  = 1200
particle = 3072**3

pi = np.pi
G  = 6.67408e-11	
H  = 100*1*1000/(3.08568E22)
rho_m = 3.0*H*H/(8.0*pi*G)

rh  = rho_m/(1.9891E30)*(3.08568E22**3)*omega_m       #(M_sun/h)/(Mpc/h)^3
print(rho_m)
print(rh)
res = rh/(particle/boxsize**3)
print(res)
