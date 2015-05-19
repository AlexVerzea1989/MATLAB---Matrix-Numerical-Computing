function [ t1,t2,t3,t4,cb] = test1 (dim)
k=4;
t1=zeros(size(dim));
t2=zeros(size(dim));
t3=zeros(size(dim));
t4=zeros(size(dim));
cb=zeros(size(dim));
j=0;
for n=dim
	j++;
	i=1:n;
	D=10.^((i-1)*(k/(n-1)));
	%D=ones(1,n);
	%for i=1:n/2
	%	D(i)=10^k;
	%end
	%D=D/(10^(k/2));
%for t=1:100
	%[Q,~]=qr(randn(n));
	%[P,~]=qr(randn(n));
	%A=Q*diag(D)*P;
	A=randn(n)*diag(D);
	B=A;
	Z=eye(n);
	condB=cond(B);
	tic;
while true
	[R,P]=QRZ_PHH(B,n);
	U=diag(diag(R))\R;
	Z1=round(inv(U));
	B=R*Z1;
	newcondB=cond(B);
	%disp([condB,newcondB]);
	if(newcondB>=condB-0.5)
		break;
	end
	Z=Z*P;
	Z=Z*Z1;
	condB=newcondB;
end
	B=A*Z;
	t=toc;
	cb(j)=cond(B);
%	disp([cond(B),cond(A)])
%	if(cond(B)<cond(A)/2)
%		break;
%	end
%end
	tic
	PLLL(A,n);
	t1(j)=toc;
	tic
	PLLL(B,n);
	t2(j)=toc+t;
	%RB=QRZ_PHH(B,n);
	tic
	LLL(A,n);
	t3(j)=toc;
	tic
	LLL(B,n);
	t4(j)=toc+t;
end;
disp([t1;t2;t3;t4;cb]);
figure;
hold on;
plot(dim,t1,'b');
plot(dim,t2,'r');
plot(dim,t3,'b:');
plot(dim,t4,'r:');
figure;
plot(dim,log10(cb));
endfunction
