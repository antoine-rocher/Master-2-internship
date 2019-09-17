import numpy as np
import subprocess

#Optimation du biais par brute force 
def compute_model(f1,f2):
    run_gs =  subprocess.Popen("/home/arocher/Stage/gaussian_stream_model/code_GSRSD/GSRSD.exe %f %f"%(f1,f2), shell=True, stdout=subprocess.PIPE).stdout
    x = np.array(run_gs.read().decode().split(" "))
    x = np.delete(x,-1)           #enleve l'espace a la fin
    data = x.astype(np.float)     #lis et convertit les outputs du GS
    r = data[data%2.5==0]
    xi0 = empty(len(r))
    xi2 = empty(len(r))
    g=1
    d=2
    for i in np.range(len(r)):
        xi0[i] = data[g]
        xi2[i] = data[d]
        g += 3
        d += 3
    return xi0[np.range(4,20)],xi2[np.range(4,20)],r[np.range(4,20)]

#Calcul du xi^2 ou x est le monopole de la simu
def Xi_2 (bias,x):
    f1,f2 = bias 
    y_model = compute_model(f1,f2)[0]
    print(f1,f2)
    return sum(y_model-x)**2
    
