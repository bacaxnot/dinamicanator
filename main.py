#%%
#EXTERNAL LIBRARIES

#INTERNAL LIBRARIES
from dynamicMethods import *;
#%%
#FINDING GEOMETRIC PROPERTIES
n = 5;
h = 3200; #m
m = 272330.6418; #kg
E = 4700*np.sqrt(27); #MPa

I1 = 300*(450**3)/12; #mm4
I2 = 450*(300**3)/12; #mm4

K1 = (12*E*I1/(h**3))*10**3; #N/m
K2 = (12*E*I2/(h**3))*10**3; #N/m

KT1 = 18*K1 + 4*K2; #horizontal
KT2 = 18*K2 + 4*K1; #vertical
#%%
#MODEL DATA
d1 = oneDeg(0, mass=m, stiffness=KT1, heigth=h)
d2 = oneDeg(1, mass=m, stiffness=KT1, heigth=h)
d3 = oneDeg(2, mass=m, stiffness=KT1, heigth=h)
d4 = oneDeg(3, mass=m, stiffness=KT1, heigth=h)
d5 = oneDeg(4, mass=m, stiffness=KT1, heigth=h)

#%%
#PROCEDURE
model = dynamicModel()
model.addDeg(d1, d2, d3, d4, d5)


w2, phi = model.modes;
wFun = np.sqrt(np.min(w2));
Tfun = 2*np.pi/wFun;
print(Tfun)


#Fundamental periods
Tx = 0.733; #s
Ty = 0.733; #s

#K determination (for EHF)
T = Tfun;
if T <= 0.5:
    ki = 1;
elif 0.5<T and T<=2.5:
    ki = 0.75 + 0.5*T;
elif T > 2.5:
    ki = 2;

spectre = aDesignSpectre(Aa=0.25,Av=0.25,Fa=1.15,Fv=1.55,I=1)
# spectre.print()

#Drifts X direction
Sax = spectre.getSa(T)
dx, vix, vx = EHF(model, Sax, ki)
print("\nHEF X DIRECTION:\n")
print(f"Base shear [kN]: {vx/1000}")
for i in range(n):
    print(f"Story {n-i}:")
    print(f"Shear [kN]: {vix[i,0]/1000}")
    print(f"Drift [%]: {dx[i,0]*100}")

#%%
