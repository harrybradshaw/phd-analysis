import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters, fit_report
import os, sys

Ic_idx = 1

def Reader(filen):
    with open(filen, 'r') as rfile:
        data = np.genfromtxt(rfile,delimiter=",")  
    #print data
    return data

def JJ_Curve(B,shift,L,d,Ic,I0):
    """
    B: mT
    shift: mT
    L: nm
    d: nm
    Ic: mA
    """
    phi_0 = 2067.8
    Lam = 45
    B_T = (B-shift) * 1E-3
    phi = B_T * L*(d + (2*Lam))
    arg = (np.pi*phi)/phi_0
    JJ =  np.abs(np.sin(arg)/arg)*Ic + I0
    return JJ

def residual(pars,x,data = None, eps=None):
    parvals = pars.valuesdict()
    shift = parvals['shift']
    L = parvals['L']
    d = parvals['d']
    Ic = parvals['Ic']
    Lam = parvals['Lam']
    I0 = parvals['I0']
    phi_0 = 2067.8
    B_T = (x-shift) * 1E-3
    phi = B_T * L*(d + (2*Lam))
    arg = (phi)/phi_0
    model = np.abs(np.sinc(arg))*Ic + I0
    #model = np.abs(np.cos(arg))*Ic + 0.28

    if data is None:
        return model
    if eps is None:
        return model-data

p_start = Parameters()
p_start.add('shift',value=0,min=-10,max=10,vary=True)
p_start.add('L',value=400,min=50,max=1000,vary=True)
p_start.add('d',value=30,min=0,vary=False)
p_start.add('Ic',value=0.02,min=0,vary=True)
p_start.add('I0',value=0,min=0.0,vary=False)
p_start.add('Lam',value=80,min=0.1,vary=False)


foldern = "/Users/harrybradshaw/OneDrive - University of Cambridge/Projects/AEN/AEN039C/JJ4 - 20210815/JJ4/IcH_T=4p2K"
filen = os.path.join(foldern,'JJ_IcH_output.txt')

#
direction = 'vir'
f_out = Reader(filen)
full_x = f_out[:,0]
full_y = f_out[:,Ic_idx]

if direction == 'neg':
    x_min = np.abs(full_x + 30).argmin()
    x_max = np.abs(full_x - 50).argmin()
elif direction == 'vir':
    x_min = np.abs(full_x + 50).argmin()
    x_max = np.abs(full_x - 50).argmin()
else:
    x_min = np.abs(full_x + 50).argmin()
    x_max = np.abs(full_x - 30).argmin()


x = f_out[x_min:x_max,0]
data = f_out[x_min:x_max,Ic_idx]
plt.plot(full_x,full_y,label='data',marker='o',markersize='3')
plt.xlim(-100,100)
B = np.linspace(-150,150,num=5001)
#popt,covar = curve_fit(JJ_Curve,data[:,0],data[:,1],p0=[11,509,35,1.65])
#JJ = JJ_Curve(B,*popt)
out = minimize(residual,p_start,args=(x,),kws={'data': data})
plt.plot(B,residual(out.params,B))
print(fit_report(out))
#print('B_min = ' + str(1/()))
a = ''
for name, par in out.params.items():
     a+=("  %s: value=%f +/- %f \n" % (name, par.value, par.stderr))

oneflux = 2067.8/(out.params["L"].value*(out.params["d"].value+(2*out.params["Lam"].value)))*1e3
print(oneflux)
#plt.plot(B9,JJ,label='sim')
direct = os.path.dirname(filen)
outname_file = os.path.join(direct,'fitting_output.txt')
outname_fig = os.path.join(direct,'fitting_output.png')
plt.ylabel(r'$I_c$ (mA)')
plt.xlabel(r'$\mu_0 H_{App}$ (mT)')
plt.ylim(0)
#plt.text(20,3,a)
plt.title(r'$I_cR_n$ = '+f'{out.params["Ic"].value * f_out[-1,4]:.5f}'+' mV')
plt.savefig(outname_fig,dpi=400)
plt.show()


np.savetxt(outname_file,np.stack((B,residual(out.params,B)),axis=1),delimiter=',')




