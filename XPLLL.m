function [Ay,Z,l,u,inactive]=XPLLL(Ay,l,u,n,inactive,Z)
%REDUCT: LLL reduction
%    LLL(A,y) Reduce ILS: min ||y-Ax||_2 using default config.
%    A: m by n matrix to be reduced
%    d2,[1,2) (1)   :sigma
%    Fqrz: initial qr strategy
%    Fpermu: mutation updating strategy
%
% Output arguments:
%    R: Out - n by n factor matrix (strict upper triangular)
%    Z: Out - n by n Unimodular matrix
%
%
%global Rho_2;
[~,ny]=size(Ay);
if nargin<4
    n=ny-1;
end

%set_count(); %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%set_count(0,2); %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% [Ay,Z]=QRZ_PHH(Ay,n);
% l=Z'*l;
% u=Z'*u;
%

% [r2,Ay,l,u,p,Z,Ry,Z2,pp]=PLLL_residual(Ay,l,u);
% R=Ay(1:n,1:n);
% y=Ay(1:n,ny);
% [inactive,l,u,r2]=inactiveset2(R,y,l,u,r2);
% 
% Rho_2=r2;
% if(all(inactive))
%     Ay=Ry;
%     %l=l(pp);
%     %u=u(pp);
%     Z=Z2(p,:);
%     return;
% end
    
if nargin<6
[R,Z,y,l,u] = QRZ_BC(Ay(:,1:n),Ay(:,ny),l,u);
Ay=[R,y];
[inactive,l,u]=inactiveset2(R,y,l,u);
end
%disp([2,sum(inactive),length(inactive)]);
%

%invZ=Z';%extra
k=2;
maxk=n;
%disp(inactive');
while k<=n
    k1=k-1;%1
    if(inactive(k1))
        mu = round(Ay(k1,k)/Ay(k1,k1)); %2
        t = Ay(k1,k)-mu*Ay(k1,k1); %2
    else
        if k>maxk
            maxk=k;
            k=k+1;
            continue;
        end
        mu=0;
        t = Ay(k1,k);
    end
    
    % add_count(10);%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    if Ay(k1,k1)^2 - (t^2+Ay(k,k)^2) > 1e-14 %5
        if mu~=0
            Ay(1:k1,k)=Ay(1:k1,k)-mu*Ay(1:k1,k1);
            Z(:,k)=Z(:,k)-mu*Z(:,k1);
            l(k1)=-inf;
            u(k1)=inf;
            [Ay,Z,l,u]=extended_IGT(Ay,Z,k,inactive,l,u);
%             i1=-1;
%             for i=k-2:-1:1
%                 if(inactive(i))
%                     i1=i;
%                     break;
%                 end
%             end
%             if i1~=-1
%                 for i=i1-1:-1:0
%                     if(i==0 || inactive(i))
%                         i2=i+1;
%                         if i1~=i2
%                             temp=Ay(i2:i1,i1);
%                             mu=round((temp'*Ay(i2:i1,k))/(temp'*temp));
%                         else
%                             mu=round(Ay(i1,k)/Ay(i1,i1));
%                         end
%                         if mu~=0
%                             Ay(1:i1,k)=Ay(1:i1,k)-mu*Ay(1:i1,i1);
%                             Z(:,k)=Z(:,k)-mu*Z(:,i1);
%                         end
%                         i1=i;
%                     end
%                 end
%             end
        end
        Ay(:,[k1,k])=Ay(:,[k,k1]);
        l([k1,k])=l([k,k1]);
        u([k1,k])=u([k,k1]);
        inactive([k1,k])=inactive([k,k1]);
        Z(:,[k1,k])=Z(:,[k,k1]);
        %invZ([k1,k],:)=invZ([k,k1],:);%extra
        
        r=sqrt(Ay(k1,k1)^2+Ay(k,k1)^2);
        c=Ay(k1,k1)/r;
        s=Ay(k,k1)/r;
        Ay(k1,k1)=r;
        Ay(k,k1)=0;
        Ay([k1,k],k:ny)=[c,s; -s,c]*Ay([k1,k],k:ny);
        %[A,Z,invZ]=Fpermu(A,Z,k,invZ);
        % add_count(1,2); %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        if k>2
            k=k-1;
        end
    else
        if k>maxk
            maxk=k;
        end
        k=k+1;
    end
end
%  [R,y,p,l,u] = QRP_NBC(Ay(:,1:n),Ay(:,ny),l,u);
%  Ay=[R,y];
%  Z=Z(:,p);
%  inactive=inactive(p);
end

function [R,Z,l,u]=extended_IGT(R,Z,k,inactive,l,u)
indx=[0;find(inactive(1:k-2))];
len=length(indx);
for i=len:-1:2
    i1=indx(i);
    i2=indx(i-1)+1;
    if i1~=i2
        temp=R(i2:i1,i1);
        mu=round((temp'*R(i2:i1,k))/(temp'*temp));
    else
        mu=round(R(i1,k)/R(i1,i1));
    end
    if mu~=0
        R(1:i1,k)=R(1:i1,k)-mu*R(1:i1,i1);
        %F(i1,k:end)=F(i1,k:end)+mu*F(k,k:end);
        Z(:,k)=Z(:,k)-mu*Z(:,i1);
        l(i1)=-inf;
        u(i1)=inf;
    end
end

end
