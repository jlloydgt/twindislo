import numpy as np
from math import pi

def sigxxstat(xx, zz, burg, matprops):
    lam, mu, rho, nu = matprops.lam, matprops.mu, matprops.rho, matprops.v
    ans = -burg*mu/(2.0*pi*(1.0-nu))*zz*(3.0*xx**2+zz**2)/(xx**2+zz**2)**2
    return np.where(np.isnan(ans), 0.0, ans)
    
def sigzzstat(xx, zz, burg, matprops):
    lam, mu, rho, nu = matprops.lam, matprops.mu, matprops.rho, matprops.v
    ans = burg*mu/(2.0*pi*(1.0-nu))*zz*(xx**2-zz**2)/(xx**2+zz**2)**2
    return np.where(np.isnan(ans), 0.0, ans)    

def sigxzstat(xx, zz, burg, matprops):
    lam, mu, rho, nu = matprops.lam, matprops.mu, matprops.rho, matprops.v
    ans = burg*mu/(2.0*pi*(1.0-nu))*xx*(xx**2-zz**2)/(xx**2+zz**2)**2
    return np.where(np.isnan(ans), 0.0, ans)

def sigstat(x, z, w, burg, matprops):
    s, c = np.sin(w), np.cos(w)

    #solve in prime coordinates
    xp =  x*c+z*s
    zp = -x*s+z*c

    sigxxp = sigxxstat(xp,zp,burg,matprops)
    sigzzp = sigzzstat(xp,zp,burg,matprops)
    sigxzp = sigxzstat(xp,zp,burg,matprops)

    #rotate stress from prime to reference - o = rt.o'.r
    sigxx = c*(sigxxp*c-sigxzp*s)-s*(sigxzp*c-sigzzp*s)
    sigzz = s*(sigxzp*c+sigxxp*s)+c*(sigzzp*c+sigxzp*s)
    sigxz = s*(sigxxp*c-sigxzp*s)+c*(sigxzp*c-sigzzp*s)    

    return sigxx, sigzz, sigxz

def sigxztiltGB(x, z, x0, angle, sign, matprops, burg=7.621554e-4, small=1e-8, threshold=5.0):
    #pp. 733 Hirthe and Lothe
    #pp. 94 Cottrell, 1953
    if angle < small:
        return np.zeros((np.shape(x))), np.zeros((np.shape(x))), np.zeros((np.shape(x)))
    else:
        lam, mu, rho, nu = matprops.lam, matprops.mu, matprops.rho, matprops.v    
        D = burg/(2.0*np.sin(angle/2.))
        X, Z = (x-x0)/D, (z%D)/D #sign*z/D
        
        #don't use np.where because it computes the value then overwrites it, which still causes an error
        ge, le = np.where(abs(X) > threshold), np.where(abs(X) <= threshold)
        sigxz = np.zeros(x.shape)
        sigxz[ge] = 2.*np.pi*mu*burg*X[ge]/(D*(1.-nu))*np.exp(-2.*np.pi*abs(X[ge]))*np.cos(2.*np.pi*Z[ge])
        sigxz[le] = mu*burg/(2.*D*(1.-nu)*(np.cosh(2.*pi*X[le])-np.cos(2.*pi*Z[le]))**2)*(2.*pi*X[le])*(np.cosh(2.*pi*X[le])*np.cos(2.*pi*Z[le])-1.)

        return sign*sigxz

'''
#NOT WORKING YET (COSH OF LARGE VALUES GIVES ERRORS)
def sigtiltGB(x, z, x0, angle, sign, matprops, burg=7.621554e-4, small=1e-8):
    #pp. 733 Hirthe and Lothe
    #pp. 94 Cottrell, 1953
    if angle < small:
        return np.zeros((np.shape(x))), np.zeros((np.shape(x))), np.zeros((np.shape(x)))
    else:
        lam, mu, rho, nu = matprops.lam, matprops.mu, matprops.rho, matprops.v    
        D = burg/(2.0*np.sin(angle/2.))
        print D
        X, Z = (x-x0)/D, sign*z/D
        np.seterr(all='print')
        try:
            sig0 = mu*burg/(2.*D*(1.-nu)*(np.cosh(2.*pi*X)-np.cos(2.*pi*Z))**2)
        except:
            print X    
        sigxx = -sig0*np.sin(2.*pi*Z)*(np.cosh(2.*pi*X)-np.cos(2.*pi*Z)+2.*pi*X*np.sinh(2.*pi*X))
        sigxz = sig0*(2.*pi*X)*(np.cosh(2.*pi*X)*np.cos(2.*pi*Z)-1.)
        sigzz = -sig0*np.sin(2.*pi*Z)*(np.cosh(2.*pi*X)-np.cos(2.*pi*Z)-2.*pi*X*np.sinh(2.*pi*X))
        
        siggxx = np.where(np.isnan(sigxx), 0.0, sigxx)
        siggxz = np.where(np.isnan(sigxz), 0.0, sigxz)
        siggzz = np.where(np.isnan(sigzz), 0.0, sigzz)
        
        return sigxx, sigzz, sigxz
'''


