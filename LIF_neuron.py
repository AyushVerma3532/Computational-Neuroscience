import math
import numpy as np
import matplotlib.pyplot as plt

class InputError(Exception):
    pass

def tospikeval(spikes, n):# converts spikes into stream of 0s and 1s
    spikeval=np.zeros(n, dtype=int)
    j=0
    if len(spikes)!=0:
        for k in range(0, n):
            if spikes[j]==k:
                j+=1
                spikeval[k]=1
            if j>=len(spikes):
                break
    return spikeval

def rate_encode(dec_val, mean_val, std_dev, density, zero_state):# encodes devcimal values between 0 and 1 to 
    out=[]
    for v in dec_val:
        if v> 1 or v<0:
          raise InputError("dec_val can contain numbers in range [0,1] only")

        rate=[]
        random_normal=np.random.normal(mean_val, std_dev, int(density*v+zero_state))
        for r in random_normal:
            if r>=0 and r<2*mean_val:
                rate.append(int(r))
        rate=np.array(rate)
        rate=np.unique(rate)
        out.append(rate)
    return out

def LIF(inputspikes, weights, tau_m, tau, vt=1, dvt=0.2, reset=0, precision=10):
    n=0
    for s in inputspikes:
        if len(s)!=0:
            if s[-1]>n:
                n=s[-1]
    final=np.zeros(n)
    for k in range(0, len(inputspikes)):
        final+=weights[k]*tospikeval(inputspikes[k], n)
    vm=0.0
    vth=vt
    j=0
    t_vm=0
    t_vt=0
    vmarr=[]
    vtarr=[]
    tarr=[]
    vm0=0
    vt0=vt
    outarr=[]
    for i1 in range(0, precision*n):
        i=i1/precision
        tarr.append(i)
        if i1%precision==0:
            if final[int(i)]!=0.0:
                t_vm=0
                vm=vm+final[int(i)]
                vm0=vm
        vm=vm0*math.exp(-t_vm*tau_m/(1000*precision))
        vm=round(vm, 3)
        vmarr.append(vm)
        if vm>=vth:
            t_vt=0
            outarr.append(int(i))
            vm0=reset
            vth+=dvt
            vt0=vth
        vth=vt+(vt0-vt)*math.exp(-t_vt*tau/(1000*precision))
        vt=round(vt, 3)
        vtarr.append(vth)
        t_vm+=0.1
        t_vm=round(t_vm, 1)
        t_vt+=0.1
        t_vt=round(t_vt, 1)
    plt.plot(tarr, vmarr, label="Membrane Potential")
    plt.plot(tarr, vtarr, label="Threshold Potential")
    plt.xlabel("time(us)")
    plt.ylabel("Potential(V)")
    plt.legend()
    plt.show()
    return np.array(outarr)
  
val=rate_encode(dec_val=[0.2, 0.2, 0.6], mean_val=50, std_dev=25, density=40, zero_state=5)
LIF(inputspikes=val, weights=[0.4, 0.2, 0.6], vt=1.0, dvt=0.2, reset=0, tau_m=8000, tau=400, precision=100)
