syms Vcc Vth Vbe1 Vbe2 Re R Rs Rb Rb1 Rb2 RL Ic1 Ic2 Ie1 Ie2 hfe1 hfe2 hfe22 gm1 gm2 hie1 hie2 hie22  rce1 rce2 GI GIS GV Ri Ria Ro R1 R2 Roa    

HFE1=501.7;                        
ALPHA1=(HFE1+1)/HFE1;
HFE2=529.8;
hfe1=HFE1;
hfe2=HFE2;
Vcc=12;
Vbe1=0.7;
Vbe2=0.75;
Rs=50;

R=33e3;
RL=15;
Re=1e3;
Rb1=470e3;
Rb2=820e3;


Vth = Vcc*Rb2/(Rb1+Rb2);
Rb = parallel(Rb1,Rb2);
Ic2=(Vth-Vbe1-Vbe2*(1/(HFE1*R)+ALPHA1*Re/R))/(Re*(1+ALPHA1/HFE2)+Rb/(HFE1*HFE2));
Ic1=Vbe2/R+Ic2/HFE2;
Vce2=Ic1*ALPHA1*Re+Ic2*Re-Vcc;
Vce1=Vcc-Ic1*(R+ALPHA1*Re)-Ic2*(-R/HFE2+Re);

Rd=parallel(Re,RL);
gm1=Ic1/25e-3;
gm2=Ic2/25e-3;
hie1=(hfe1+1)/gm1;
hie2=(hfe2+1)/gm2;
rce1=52.64/Ic1;
rce2=26.03/Ic2;

R1=parallel(Rs,Rb)+hie1;
R2=parallel(hie2,R);
hfe22=hfe2*R/(R+hie2);

GI=(hfe1*hfe2*R/(hie2+R)+hfe1)*Re/(RL+Re)
Ri=hie1+GI*Rd;
Ria=parallel(Ri,Rb)
GV=GI*Rd/(GI*Rd+hie1)
GIS=GI*Rb/(Ri+Rb)

Ro=1/(1/R1+(hfe22+1)*(rce1*hfe1/R1-1)/(rce1+R2));
Roa=parallel(Ro,rce2)


function para = parallel(x,y)
    para = (x*y)/(x+y); 
end