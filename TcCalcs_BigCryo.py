
"""
Created on Wed Nov 14 21:21:11 2018

@author: harrybradshaw
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

def reader(filen):
    with open(filen, 'r') as rfile:
        headers = (rfile.readline()).rstrip('\r\n').split('\t')
        try:
            r_idx = headers.index('R_nv')
        except ValueError:
            r_idx = 2

        try:
            B = headers.index('B_digital_(T)')
        except:
            B = 1
        try:
            T_idx = headers.index('sensor_B_(K)')
        except ValueError:
            T_idx = 1
        data = np.genfromtxt(rfile,skip_header=0)

    return data, r_idx, T_idx, B

def TcPlot(filen):
    data, r_idx, T_idx, B = reader(filen)
    idx_10 = np.nanargmin(np.abs(data[:,T_idx]- 9))
    if idx_10 > 2:
        r_10 = np.mean(data[(idx_10-2):(idx_10 +2),r_idx])
    else:
        r_10 = np.mean(data[(idx_10):(idx_10 +2),r_idx])


    # r_10 = 4.52089
    ndata = data
    ndata[:,r_idx] = data[:,r_idx] / r_10

    plt.plot(ndata[:,T_idx],ndata[:,r_idx])
    plt.xlabel('T (K)', fontsize=20)
    plt.ylabel('Normalised Resistance', fontsize=20)
    # plt.legend()
    up_b = 0.9
    low_b = 0.1

    idx_high = np.nanargmin(np.abs(ndata[:,r_idx]-up_b))
    idx_low = np.nanargmin(np.abs(ndata[:,r_idx]-low_b))
    idx_mid = np.nanargmin(np.abs(ndata[:,r_idx]-(0.5)))

    print(r_10)

    tc_width = ndata[idx_high,T_idx] - ndata[idx_low,T_idx]
    tc = ndata[idx_high,T_idx]
    plt.plot(tc,ndata[idx_high,r_idx],marker='o')
    print(ndata[idx_high,T_idx])
    plt.plot(ndata[idx_mid,T_idx],ndata[idx_mid,r_idx],marker='o')
    print(ndata[idx_mid,T_idx],tc_width)
    plt.plot(ndata[idx_low,T_idx],ndata[idx_low,r_idx],marker='o')
    print(ndata[idx_low,T_idx])
    saven = filen.split('.')[0] + '.png'
    plt.tight_layout()
    plt.savefig(saven,dpi=400)

def linear(x,m,c):
    return (m*x)+c

def norm_funct(x,min_x,max_x,a=-1,b=1):
    """
    Normalises data in the interval [a,b]
    """
    return (b-a)*((x-min_x)/(max_x - min_x)) + a

def TcLinear(filen, T_ref = 10):
    print(filen)
    data, r_idx, T_idx, B = reader(filen)
    #print(f'B_av = {np.mean(data[:,B])}')
    idx_ref = np.nanargmin(np.abs(data[:, T_idx] - 9))

    if idx_ref > 2:
        r_ref = np.mean(data[(idx_ref - 2):(idx_ref + 2), r_idx])
    else:
        r_ref = np.mean(data[idx_ref:(idx_ref + 2), r_idx])

    #r_ref = 0.1652055
    data[:, r_idx] = norm_funct(data[:, r_idx], 0, r_ref)
    print(r_ref)
    plt.plot(data[:,r_idx],data[:,T_idx],marker='o',markersize=1)

    up_b = 0.5
    low_b = -0.5
    mid_b = 0.0

    idx_high = np.nanargmin(np.abs(data[:, r_idx] - up_b))
    idx_low = np.nanargmin(np.abs(data[:, r_idx] - low_b)) + 1
    idx_mid = np.nanargmin(np.abs(data[:, r_idx] - mid_b))

    #idx_low = 201
    if(idx_low-idx_high == 1):
        idx_low = idx_high \
                  + 2
    print(idx_high,idx_low)


    popt, covar = curve_fit(linear, data[idx_high:idx_low, r_idx], data[idx_high:idx_low, T_idx])
    err = np.sqrt(np.diag(covar))
    T_vals = np.linspace(-0.5,0.5,num=10)
    plt.plot(T_vals,linear(T_vals, *popt))
    plt.plot( data[idx_high:idx_low, r_idx], data[idx_high:idx_low, T_idx],marker='o',markersize=1)
    print(f'Tc = {popt[1]}\t{err[1]}\t{popt[0]}')
    print(f'Gradient = {popt[0]}')
    print(f'Extrap Tc_width = {popt[0]*1.6}')
    #plt.ylim(6.8,9.2)
    #plt.xlim(-1,1)


file1 = "/Users/harrybradshaw/OneDrive - University of Cambridge/Projects/NHJ/NHJ052/52A/Data/OOP_Transport/"
file2 = "52A_OOP_Tc_Cool_dir=neg_set=-0p80T_cool.txt"

file1 = "/Users/harrybradshaw/OneDrive - University of Cambridge/Projects/NHJ/NHJ063/IP_Transport"
file2 = "63A_dir=neg_set=-0p50T_cool.txt.txt"

def filenamer(file1,file2):
    filen = os.path.join(file1, file2)
    return filen

def looper():
    for i in range(8):
        st = "set="+str(i)+"p0_dir=vir.dat"
        try:
            print(st)
            TcLinear(filenamer(file1,st))
        except:
            print(st+" not found")



#looper()
TcLinear(filenamer(file1,file2))
plt.show()