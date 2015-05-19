n=50;
m=60;
[A,y,x]=randILS(m,n);
PRy=PLLL([A,y],n,1,@QRZ_HH,@permu_Givens);
Ry=LLL([A,y],n,1,@QRZ_HH,@permu_Givens);
d1=diag(PRy);
d2=diag(Ry);
max(abs(d1-d2))
