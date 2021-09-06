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
    return np.piecewise(I,[I<=Ic,I>Ic],[lambda I: g, lambda I: Ic*Rn*np.sqrt((I/Ic)**2 - 1) + g])

def IcFit_neg(I,Ic,Rn,g):
    #return Ic*Rn*np.sqrt((I/Ic)**2 - 1) * (I<Ic)
    return np.piecewise(I,[I>=-Ic,I<-Ic],[lambda I: g, lambda I: Ic*Rn*np.sqrt((I/Ic)**2 - 1) + g])

        
filen = "/Users/harrybradshaw/OneDrive - University of Cambridge/Projects/NHJ/NHJ042/Data/20201221 Harry 42B/JJ3/IcH_pos_to_neg_T=4p0K_set=1p0T/0.150.dat"

def Fitter(filen,direct,plot_indiv=False):
    field = float(((filen.split(direct)[1]).rstrip('.dat')).lstrip('/'))
    data = Reader(filen)
    #data[:,2] += (1.43783249e-05)
    upper_thres = 3
    dVdI_thres = 1
    count = 0
    points = len(data[:,0])
    for i in range(points):
        if data[i,0] < upper_thres:
            count +=1
        if count > 1:
            break
    data = data[i-count+2:points-i+count-2,:]
    points = len(data[:,0])
    half = int(points/2.0)

    data_pos = data[half:,:]
    data_neg = data[:half,:]

    count = 0
    for i in range(len(data_neg)):
        if data_neg[i,0] < dVdI_thres:
            count=1
            break
    if(count):
        Ic_neg_guess = np.abs(data_neg[i,1])
    else:
        Ic_neg_guess = 0

    for i in range(len(data_pos)-1,-1,-1):
        if data_pos[i,0] < dVdI_thres:
            break
    Ic_pos_guess = np.abs(data_pos[i,1])
   
    print(Ic_neg_guess,Ic_pos_guess)
    initial = [Ic_neg_guess,2.1,0]
    initial_pos = [Ic_pos_guess,2.1,0]

    try:
        popt_p,covar = curve_fit(IcFit_pos,data_pos[:,1],data_pos[:,2],p0=initial_pos)
        popt_n,covar = curve_fit(IcFit_neg,data_neg[:,1],data_neg[:,2],p0=initial)
    except RuntimeError:
        popt_p = [0,0,1]
        popt_n = [0,0,1]

    if plot_indiv == True:
        plt.plot(data[:,1],data[:,2])
        plt.plot(data_pos[:,1],IcFit_pos(data_pos[:,1],*popt_p))
        plt.plot(data_neg[:,1],IcFit_neg(data_neg[:,1],*popt_n))
        plt.xlabel('I')
        plt.ylabel('V')
        plt.show()
        

    Icn = popt_n[0]
    Icp = popt_p[0]
    Rn = ((np.abs(popt_p[1]) + np.abs(popt_n[1]))/2.0)
    mid = data[half,0]
    print(field,Icn,Icp,Rn)
    return Icn,Icp,Rn,field,mid,Ic_neg_guess,Ic_pos_guess

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

def Directory(direct,debug_plot=False):
    Icn = []
    Icp = []
    IcAv = []
    Rn = []
    B = []
    Mid = []
    IcPg = []
    IcNg = []
    files = GatherFiles(direct)
    print('Temp \t Icn \t Icp \t Rn')
    for fname in files:
        try:
            icn,icp,rn,fe,mid,icng,icpg = Fitter(fname,direct,plot_indiv=debug_plot)
            
        except (ValueError, TypeError) as e:
            icn,icp,rn = 0,0,0
            fe = float(((fname.split(direct)[1]).rstrip('.dat')).lstrip('/'))
            mid = 0
        Icn.append(icn*1e3)
        Icp.append(icp*1e3)
        IcAv.append(((icn+icp)/2.0)*1e3)
        Rn.append(rn)
        B.append(fe)
        Mid.append(mid)
        IcPg.append(icpg*1e3)
        IcNg.append(icng*1e3)
    

    #plt.plot(B,np.array(Icn))
    #plt.plot(B,np.array(Icp))
    plt.plot(B,np.array(IcPg))
    plt.plot(B,np.array(IcNg))
    plt.ylabel('Ic (mA)')
    plt.ylim(0)
    plt.show()
    plt.plot(B,Mid)
    plt.show()
    print(direct)
    print(np.mean(np.array(Rn)))
    outname = os.path.join(direct,'IcT_output.txt')
    header = 'T,Ic,Ic+_Fit,Ic-_Fit,Ic+_Guess,Ic-_Guess,Rn_Fit,R(I=0)\nK,mA,mA,mA,mA,mA,Ohm,Ohm'
    np.savetxt(outname,np.stack((B,IcAv,Icp,Icn,IcPg,IcNg,Rn,Mid),axis=1),delimiter=',',header=header)

if(__name__=="__main__"):
    Directory("/Users/harrybradshaw/OneDrive - University of Cambridge/Projects/Guang_CrossJunction/IcT_new",debug_plot=False)
    #Fitter('/Users/harrybradshaw/OneDrive - University of Cambridge/Projects/AEN/AEN039A/JJ4 - Remeasure/IcT/6.477.dat','Users/harrybradshaw/OneDrive - University of Cambridge/Projects/AEN/AEN039A/JJ4 - Remeasure/IcT',plot_indiv=True)