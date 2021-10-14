from classes import *;

def EHF(model:dynamicModel,
        sa:float,
        ki:float,
        g:float = 9.81):
    """
    Equivalent Horizontal Force Method

    Arguments:
    - model: target dynamicModel object
    - sa: spectre acceleration
    """
    #Individual properties
    m = np.array([element.m for element in model.elements]);
    h = np.array([element.h for element in model.elements])/1000;
    k = model.kMatrix
    n = len(m); #Number of stories
    #EHF: global and per element
    v = sa*g*np.sum(m);
    for i in range(n):
        if i == 0:
            pass
        else:
            h[i] = h[i] + h[i-1]
    h = np.flip(h)
    vi = ( v*(m*(h**ki))/np.sum(m*(h**ki)) ).reshape((n, 1));
    #Displacements (dxi) and drifts (di)
    dxi = np.linalg.solve(k, vi);
    for i in range(n):
        if i == n-1:
            pass
        else:
            dxi[i, 0] = dxi[i, 0] - dxi[i+1, 0]
    di = dxi/h.reshape((n,1));
    return di, vi, v