function [zf,beta2] = search_fast_n(Ry,n,ny,N,delta)

[~,nn] = size(Ry);
% Check input arguments
if nargin < 2 % input error
    n=nn-1;
end
if nargin <3
    ny=n+1;
end
if(nargin<4)
    N=n;
end
if(nargin<5)
    delta=0;
end

beta2=inf;
k=n;
z=zeros(n,1);
c=zeros(n,1);
d=zeros(n,1);
s=zeros(n+1,1);
b2=zeros(n,1);
bb2=zeros(n-1,1);
S=zeros(n-1,1);
rii2=diag(Ry);
rii2=rii2(1:n).^2;
srii2=zeros(n,1);
for i=1:n-1
    srii2(i+1)=srii2(i)+rii2(i);
end
c(n)=Ry(n,ny)/Ry(n,n);
z(n)=round(c(n));
d(n)=sign(c(n)-z(n));
b2(n)=rii2(n).*(z(n)-c(n)).^2;
s(n)=s(n+1)+b2(n);
for j=1:N
    k1=k-1;
    stop=2;
    for i=k1:-1:1
        c(i)=(Ry(i,ny)-Ry(i,i+1:n)*z(i+1:n))/Ry(i,i);
        z(i)=round(c(i));
        b2(i)=rii2(i).*(z(i)-c(i)).^2;
        s(i)=s(i+1)+b2(i);
        if(s(i)>beta2)
            stop=i+1;
            break;
        end
    end
    if(stop>n)
        break;
    end
    d(stop:k1)=sign(c(stop:k1)-z(stop:k1));
    bb2(stop-1:k1)=rii2(stop:k).*(z(stop:k)+d(stop:k)-c(stop:k)).^2;
    

    
    if(s(1)<beta2)
        beta2=s(1);
        zf=z;
    end
    
    S(stop-1:k1)=(beta2-bb2(stop-1:k1)-s(stop+1:k+1))./srii2(stop:k);
    %[S(k1),beta2,bb2(k1),s(k+1),srii2(k)]
    [Sk,k]=max(S(stop-1:n-1));
    if(Sk<=delta)
        %disp([j,N]);
        break;
    end
    k=k+stop-1;
    z(k)=z(k)+d(k);
    d(k)=-d(k)-sign(d(k));
    b2(k)=rii2(k).*(z(k)-c(k)).^2;
    s(k)=s(k+1)+b2(k);
end
