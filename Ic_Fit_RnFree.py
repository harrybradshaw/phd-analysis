import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import os

g = 4*1E-5
#Rn = 0.0417

def Reader(filen):
    with open(filen, 'r') as rfile:
        data = np.genfromtxt(rfile,delimiter=",")  
    #print data
    return data

def IcFit_pos(I,Ic,Rn):
    #return Ic*Rn*np.sqrt((I/Ic)**2 - 1) * (I<Ic)
    return np.piecewise(I,[I<Ic,I>=Ic],[lambda I: g, lambda I: Ic*Rn*np.sqrt((I/Ic)**2 - 1) + g])

def IcFit_neg(I,Ic,Rn):
    #return Ic*Rn*np.sqrt((I/Ic)**2 - 1) * (I<Ic)
    return np.piecewise(I,[I>=-Ic,I<-Ic],[lambda I: g, lambda I: (-1)*Ic*Rn*np.sqrt((I/Ic)**2 - 1) + g])

        
filen = "/Users/harrybradshaw/OneDrive - University of Cambridge/JJ5/IcH_pos_to_neg_T=1p6K_set=1p0T/JJ_IcH_output.txt"
direct = '/Users/harrybradshaw/OneDrive - University of Cambridge/JJ5'
def Fitter(filen,direct=direct,debug_plot=True,generate_file=False):
    try:
        field = float(((filen.split(direct)[1]).rstrip('.dat')).lstrip('/'))
    except:
        field= 0
    data = Reader(filen)
    #data[:,2] += (1.43783249e-05)

    
    upper_thres = 0.2
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
    if(debug_plot):
        plt.plot(data[:,1],data[:,2])
    initial = [0.02E-3,0.05]
    initial_pos = [0.02E-3,0.05]
    #data_pos[:,1] = data_pos[:,1]*(-1)
    #data_neg[:,1] = data_neg[:,1]*(-1)

    global g
    g = data_pos[0,2]

    try:
        popt_p,covar = curve_fit(IcFit_pos,data_pos[:,1],data_pos[:,2],p0=initial_pos)
        if(debug_plot):
            plt.plot(data_pos[:,1],IcFit_pos(data_pos[:,1],*popt_p))
    
    except (ValueError, RuntimeError) as e:
        print('Didnt work')
        popt_p = [0]*4
    
    try:
        popt_n,covar = curve_fit(IcFit_neg,data_neg[:,1],data_neg[:,2],p0=initial)
        if(debug_plot):
            plt.plot(data_neg[:,1],IcFit_neg(data_neg[:,1],*popt_n))
    except (ValueError, RuntimeError) as e:
        print('Didnt work')
        popt_n = [0]*4

    if(debug_plot):
        plt.show()
    Icn = popt_n[0]
    Icp = popt_p[0]
    Ica = (popt_n[0] + popt_p[0])/2.0
    Rna = ((np.abs(popt_p[1]) + np.abs(popt_n[1]))/2.0)
    Rnn = np.abs(popt_n[1])
    Rnp = np.abs(popt_p[1])
    print(field,Ica*1e3,Icn*1e3,Icp*1e3,Rnn,Rnp,Rnp*Icp)

    if(generate_file):
        I = np.concatenate((data_neg[:,1],data_pos[:,1]))
        V = np.concatenate((IcFit_neg(data_neg[:,1],*popt_n),IcFit_pos(data_pos[:,1],*popt_p)))
        #p = np.stack((data_pos[:,1],IcFit_pos(data_pos[:,1],*popt_p)))
        #n = np.stack((data_neg[:,1],IcFit_neg(data_neg[:,1],*popt_n)))
        c = np.stack((I,V),axis=1)
        outname = os.path.join(direct,'IV-Fit.txt')
        np.savetxt(outname,c,delimiter=',')

    return Icn,Icp,Rnn,Rnp,Rna,field

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
    Rnn = []
    Rnp = []
    Rn = []
    B = []
    files = GatherFiles(direct)
    outname_fig = os.path.join(direct,'JJ_IcH_Figure.jpg')
    outname = os.path.join(direct,'JJ_IcH_output.txt')
    for fname in files:
        icn,icp,rnn,rnp,rna,fe = Fitter(fname,direct,debug_plot=debug_plot)
        Icn.append(icn*1e3)
        Icp.append(icp*1e3)
        IcAv.append(((icn+icp)/2.0)*1e3)
        Rnn.append(rnn)
        Rnp.append(rnp)
        Rn.append(rna)
        B.append(fe*1e3)
    
    plt.plot(B,Icn,label='Ic-')
    plt.plot(B,Icp,label='Ic+')
    plt.legend()
    plt.xlabel('Applied Field (mT)')
    plt.ylabel('Ic (mA)')
    plt.ylim(0)
    try:
        plt.savefig(outname_fig,dpi=40)
    except PermissionError:
        pass
    plt.show()
    #plt.plot(B,Rn)
    #plt.show()
    print(direct)
    print(np.mean(np.array(Rn)))
    
    np.savetxt(outname,np.stack((B,IcAv,Icp,Icn,Rn,Rnp,Rnn),axis=1),delimiter=',')

if(__name__=="__main__"):
    Directory("/Users/harrybradshaw/OneDrive - University of Cambridge/Projects/AEN/AEN039C/JJ4 - 20210815/JJ4/IcH_T=6p0K_neg_to_pos",debug_plot=False)
    #Fitter('/Users/harrybradshaw/OneDrive - University of Cambridge/Projects/AEN/AEN039A - Oven Baked/JJ4/IV_0.dat','/Users/harrybradshaw/OneDrive - University of Cambridge/Projects/AEN/AEN039A - Oven Baked/JJ4',debug_plot=True,generate_file=True)