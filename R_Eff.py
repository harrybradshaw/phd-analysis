import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import os

def Reader(filen):
    with open(filen, 'r') as rfile:
        data = np.genfromtxt(rfile,delimiter=",")  
    #print data
    return data

def IcFit_pos(I,Ic,Rn,g):
    #return Ic*Rn*np.sqrt((I/Ic)**2 - 1) * (I<Ic)
    return np.piecewise(I,[I<Ic,I>=Ic],[lambda I: g, lambda I: Ic*Rn*np.sqrt((I/Ic)**2 - 1) + g])

def IcFit_neg(I,Ic,Rn,g):
    #return Ic*Rn*np.sqrt((I/Ic)**2 - 1) * (I<Ic)
    return np.piecewise(I,[I>=-Ic,I<-Ic],[lambda I: g, lambda I: Ic*Rn*np.sqrt((I/Ic)**2 - 1) + g])

        
#filen = "/Users/harrybradshaw/OneDrive - University of Cambridge/Projects/NHJ/NHJ042/Data/20201221 Harry 42B/JJ3/IcH_pos_to_neg_T=4p0K_set=1p0T/0.150.dat"

def Fitter(filen,direct):
    field = float(((filen.split(direct)[1]).rstrip('.dat')).lstrip('/'))
    data = Reader(filen)
    #data[:,2] += (1.43783249e-05)
    I = 3E-3
    i_idx = np.argmin(np.abs(data[:,1] - I))
    print(data[i_idx,:])
    print(data[-i_idx-1,:])
    deltaI = np.abs(data[i_idx,1])+np.abs(data[-i_idx-1,1])
    deltaV = np.abs(data[i_idx,2])+np.abs(data[-i_idx-1,2])
    Ref = deltaV/deltaI
    
    return Ref,field

def GatherFiles(direct):
    files = glob.glob('{}/*.dat'.format(direct))
    p=[]
    for fname in files:
        p_ = (fname.split(direct)[1]).rstrip('.dat').lstrip('/')
        p_ = float(p_)
        p.append(p_) 

    q = sorted(p)
    r = ["{:.3f}".format(i) for i in q]
    s = [direct+'/'+str(i)+'.dat' for i in r]
    
    files = s
    
    return files

def Directory(direct):
    R = []
    B = []
    files = GatherFiles(direct)
    outname_fig = os.path.join(direct,'R_eff_Figure.jpg')
    outname = os.path.join(direct,'R_eff_output.txt')
    for fname in files:
        ref,fe = Fitter(fname,direct)
        R.append(ref)
        B.append(fe*1e3)
    
    plt.plot(B,R,label='Ic-')
    plt.legend()
    plt.xlabel('Applied Field (mT)')
    plt.ylabel('R')
    #plt.ylim(0)
    #plt.savefig(outname_fig,dpi=400)
    plt.show()
    #plt.plot(B,Rn)
    #plt.show()
    print(direct)
    #print(np.mean(np.array(Rn)))
    
    np.savetxt(outname,np.stack((B,R),axis=1),delimiter=',')

if(__name__=="__main__"):
    Directory("/Users/harrybradshaw/OneDrive - University of Cambridge/Projects/NHJ/NHJ044/Data/20210220 Harry 44B/JJ2/IcH_inc_T=4p0K_set=1p0T_90degree")