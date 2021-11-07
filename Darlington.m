syms Vcc Vth Vbe1 Vbe2 Re R Rs Rb Rb1 Rb2 RL Ic1 Ic2 Ie1 Ie2 hfe1 hfe2 hfe22 gm1 gm2 hie1 hie2 hie22  rce1 rce2 GI GIS GV Ri Ria Ro Ro1 Ro2 Roa    
HFE1=501.6;             
HFE2=488.7;
hfe1=HFE1;
hfe2=HFE2;
Vcc=12;
Vbe1=0.7;
Vbe2=0.7;
Rs=50;

R=33e3;
RL=10;
Re=1e3;
Rb1=470e3;
Rb2=820e3;

Vth = Vcc*Rb2/(Rb1+Rb2);
Rb = parallel(Rb1,Rb2);
Ie2=(Vth-Vbe1-Vbe2*(1+Re/R+Rb/(R*(HFE1+1))))/(Re+Rb/((HFE2+1)*(HFE1+1)));
Ie1=Ie2/(HFE2+1)+Vbe2/R;
Ic2=Ie2*HFE2/(HFE2+1);
Ic1=Ie1*HFE1/(HFE1+1);
Vce2=Vcc-Ie2*Re;
Vce1=Vcc-Vbe2-Ie2*Re;

Rd=parallel(Re,RL);
gm1=Ic1/25e-3;
gm2=Ic2/25e-3;
hie1=(hfe1+1)/gm1;
hie2=(hfe2+1)/gm2;
hfe22=hfe2*R/(hie2+R);
hie22= parallel (hie2,R);
rce1=52.64/Ic1;
rce2=62.64/Ic2;

GI=(hfe1+1)*(hfe22+1)*Re/(RL+Re)
Ri=hie1+(hfe1+1)*hie22+(hfe22+1)*(hfe1+1)*Rd;
Ria=parallel(Ri,Rb)
GV=(hfe1+1)*(hfe22+1)*Rd/Ri
GIS=GI*Rb/(Ri+Rb)

Ro1=(parallel(Rs,Rb)+hie1)/(1+hfe1);
Ro2=(parallel(Ro1,rce1)+hie22)/(hfe22+1);
Ro=parallel(Ro2,rce2);
Roa=parallel(Ro,Re)

function para = parallel(x,y)
    para = (x*y)/(x+y); 
end