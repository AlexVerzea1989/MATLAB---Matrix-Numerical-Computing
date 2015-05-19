function [R,by,l,u,Z]=UBILS(A,y,l,u)
if nargin==0
   [A,y,l,u,x]=randILS_box(0,[10,20],0.1,2);
end
% [m,n]=size(A);
% Ay=[A,y];
% [Ry,p]=QRP_PHH(Ay,n,@max);
% R=Ry(:,1:n);
% by=Ry(:,n+1);
% l=l(p);
% u=u(p);

%x0=max(l,min(u,A\y));
x0=(l+u)/2;
%xs=gradproj(x0,A'*A,-2*A'*y,u,l,true);
xs2=gradproj(x0,A'*A,-(A'*y),u,l,false);
%x1=max(l,min(u,round(xs)));
x2=max(l,min(u,round(xs2)));
norm(A*x-y)
%norm(A*x1-y)
norm(A*x2-y)
