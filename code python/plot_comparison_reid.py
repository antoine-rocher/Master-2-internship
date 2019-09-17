import numpy as *
import matplotlib as *

d2=loadtxt("../CLPT_GSRSD-master/data/xi_hearin_reid.txt")
v12=loadtxt("../CLPT_GSRSD-master/data/v12_hearin_reid.txt")
sp12=loadtxt("../CLPT_GSRSD-master/data/s12_hearin_reid.txt")
fxi=loadtxt("/home/arocher/Stage/gaussian_stream_model/code_GSRSD/hearin/Xi_r_CLPT_hearin.dat")
fv12=loadtxt("/home/arocher/Stage/gaussian_stream_model/code_GSRSD/hearin/V12_CLPT_hearin.dat")
fs12=loadtxt("/home/arocher/Stage/gaussian_stream_model/code_GSRSD/hearin/Sigma_CLPT_hearin.dat")

# Xi 
xtest1=fxi[:,0]
ytest1=fxi[:,1]+fxi[:,2]+fxi[:,3]+fxi[:,4]+fxi[:,5]+fxi[:,6]

#V12
xtestv12=fv12[:,0]
ytestv12=fv12[:,1]+fv12[:,2]+fv12[:,3]+fv12[:,4]+fv12[:,5]

#S12
xtests12=fs12[:,0]
ytests12par=fs12[:,1]+fs12[:,2]+fs12[:,3]+fs12[:,4]
ytests12per=fs12[:,5]+fs12[:,6]+fs12[:,7]+fs12[:,8]

#reid
xv12=v12[:,0]
yv12=v12[:,2]+v12[:,3]+v12[:,4]+v12[:,5]+v12[:,6]
xsp12=sp12[:,0]
ysp12par=sp12[:,1]+sp12[:,2]+sp12[:,3]+sp12[:,4]
ysp12per=sp12[:,5]+sp12[:,6]+sp12[:,7]+sp12[:,8]
x1=d2[:,0]
y1=d2[:,2]+d2[:,3]+d2[:,4]+d2[:,5]+d2[:,6]+d2[:,7]



figure(figsize=(50,15))
subplot(1,4,1)
loglog(xtest1,ytest1,color="red",label="my code")
loglog(x1,y1,'--',label="Wang 2013")
xlabel("Distance $r$ ($Mpc/h)$",size=20)
ylabel("$\\xi_r(r)$", size=20)
title('Correlation function in real space',size=20)
legend(fontsize=20)
subplot(1,4,2)
plot(xtestv12,ytestv12,color="red",label="my code")
plot(xv12,yv12,'--',label="Wang 2013")
xlabel("Distance $r$ ($Mpc/h)$",size=20)
ylabel("$V_{12} \ (1/f*Mpc/h)$", size=20)
title('Pairwise velocity',size=20)
legend(fontsize=20)
subplot(1,4,3)
plot(xtests12,ytests12par,color="red",label="my code")
plot(xsp12,ysp12par,'--',label="Wang 2013")
xlabel("Distance $r$ ($Mpc/h)$",size=20)
ylabel("$\sigma_{12,\parallel} \ (1/f*Mpc/h)$", size=20)
title('Velocity dispertion parallel component',size=20)
legend(fontsize=20)
subplot(1,4,4)
plot(xtests12,ytests12per*2,color="red",label="my code")
plot(xsp12,ysp12per,'--',label="Wang 2013")
xlabel("Distance $r$ ($Mpc/h)$",size=20)
ylabel("$\sigma_{12,\perp} \ (1/f*Mpc/h)$", size=20)
title('Velocity dispertion perpendicular component',size=20)
legend(fontsize=20)
#savefig('/home/arocher/Stage/gaussian_stream_model/code_GSRSD/hearin/comp_moment_CLPT_hearin.png')
