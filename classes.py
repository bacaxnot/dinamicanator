#%%
#EXTERNAL LIBRARIES
from typing import Any
from scipy import linalg as la;
import numpy as np
import matplotlib.pyplot as plt;
#INTERNAL LIBRARIES
#%%
#CLASSES
class oneDeg:
    def __init__(self,  index:int,
                        mass:float,
                        stiffness:float,
                        damping:float = 0,
                        heigth:float = 0) -> Any:
        self.index = index;
        self.m = mass;
        self.k = stiffness;
        self.c = damping;
        self.h = heigth;

class dynamicModel: #This is just for 2 FreeDeg

    def __init__(self) -> Any:
        self.elements = [];
    
    def sortElements(self):
        self.elements.sort(key=lambda x: x.index)
        
    def addDeg(self,    *newDeg:oneDeg):
        for deg in newDeg:
            self.elements.append(deg)
        self.sortElements()
    
    def setInConditions(self,   u:list = [[0],[0]],
                                v:list = [[0],[0]]):
        self.uo = np.array(u);
        self.vo = np.array(v);
    
    def setLoads(self,  p:list):
        #either [0] for free vib. or [Po, omega] for armonic forced vib.
        self.p = np.array(p);

    @property
    def kMatrix(self):
        n = len(self.elements);
        k = np.zeros((n, n));

        for i in range(n):
            if i != n-1:
                k[i:i+2, i:i+2] += self.elements[i].k*np.array([[1, -1],[-1, 1]])
            else:
                k[i,i] += self.elements[i].k

        self.k = k;
        return self.k
    @property
    def mMatrix(self):
        m = [element.m for element in self.elements]
        self.m = np.diag(m);
        return self.m
    @property
    def modes(self):
        m = self.mMatrix;
        k = self.kMatrix;

        #Getting w2 and phi
        w2, phi = la.eig(k, m);

        #Checking w2 is not complex
        isW2Complex = all(np.iscomplex(w2i) for w2i in w2);
        if not isW2Complex:
            w2 = w2.astype(float)
        else:
            print("Error, las frecuencias son números complejos")
            exit()

        #Scaling phi values
        for i in range(phi.shape[1]):
            phi[:,i] /= np.min(np.abs(phi[:,i]));

        return w2, phi

    def disengage(self):
        w2, phi = self.modes();

        self.m = phi.T @ self.m @ phi;
        self.k = phi.T @ self.k @ phi;

        return self.m, self.k, phi
#BELOW THIS COMMENT IT'S JUST EXPERIMENTAL
    def getDisplacementFunctions(self):
        m, k, phi = self.disengage()
        n = phi.shape[1]
        u = [];
        #Individual responses
        for i in range(m.shape[0]):
            #Variables
            wi = np.sqrt(k[i,i]/m[i,i])
            Ai = self.uo[i,0];
            Po = self.p[i,0];

            if Po:
                #Armonic forced vibration
                omega = self.p[i,1];
                betha = omega/wi;
                coef = Po/(k[i,i]*(1 - betha**2));
                Bi = (-coef*omega + self.vo[i,0])/wi;
                fun = lambda t: Ai*np.cos(wi*t) + Bi*np.sin(wi*t) + coef*np.sin(omega*t);
            else:
                #Free vibration
                Bi = self.vo[i,0]/wi;
                fun = lambda t: Ai*np.cos(wi*t) + Bi*np.sin(wi*t); 
            u.append(fun)
        #Modalized responses
        uc = tuple(u)
        um = [];
        um.append(lambda t : phi[0,0]*uc[0](t) + phi[0,1]*uc[1](t))
        um.append(lambda t : phi[1,0]*uc[0](t) + phi[1,1]*uc[1](t))
        return um

    def drawDisplacement(self, interval: np.arange):
        plt.figure()
        u = self.getDisplacementFunctions();
        n = len(u)
        
        for i in range(n):
            plt.subplot(n,1,i+1)
            plt.plot(interval, u[i](interval))
            plt.title(f"Masa {i+1}")
            plt.xlabel("t")
            plt.ylabel("u(t)")

        plt.show()

class aDesignSpectre:
    def __init__(self,  Aa:float,
                        Av:float,
                        Fa:float,
                        Fv:float,
                        I:float) -> Any:
        self.TO = 0.1*Av*Fv/(Aa*Fa);
        self.TC = 0.48*Av*Fv/(Aa*Fa);
        self.TL = 2.4*Fv;
        self.Aa = Aa;
        self.Av = Av;
        self.Fa = Fa;
        self.Fv = Fv;
        self.I = I;

    def firstSegment(self):
        n = round((self.TC - 0)/0.1)
        t = np.linspace(0, self.TC, n);
        Sa = 2.5*self.Aa*self.Fa*self.I*np.ones_like(t);
        return t, Sa

    def secondSegment(self):
        n = round((self.TL - self.TC)/0.01)
        t = np.linspace(self.TC, self.TL, n);
        Sa = lambda T: 1.2*self.Av*self.Fv*self.I/T;
        return t, Sa(t)

    def thirdSegment(self):
        n = round((1.5)/0.01)
        t = np.linspace(self.TL, self.TL + 1.5, n);
        Sa = lambda T: 1.2*self.Av*self.Fv*self.I*self.TL/(T**2);
        return t, Sa(t)
    
    def print(self):
        fig, ax = plt.subplots()
        ax.set_ylabel('$S_a$\n(g)')
        ax.set_xlabel('T (s)')
        ax.set_title('Acceleration Design Spectre')
        
        ax.plot(*self.firstSegment())
        ax.plot(*self.secondSegment())
        ax.plot(*self.thirdSegment())
        plt.show()

    def getSa(self, T):
        if T < 0:
            print("Error, period must be positve number")
        elif T <= self.TC:
            return 2.5*self.Aa*self.Fa
        elif T <= self.TL:
            return 1.2*self.Av*self.Fv*self.I/T
        else:
            return 1.2*self.Av*self.Fv*self.I*self.TL/(T**2)
# %%
