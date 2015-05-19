function [A,Z,l,u] = XHReduce(A,n,l,u)
%first version, slow
%  [R,y,Z] = otherConstrainedReduction(A(1:n,1:n),A(1:n,n+1));
%  A=[R,y];
[A,Z]=QRZ_HH(A,n);
for i=n:-1:1
    [a,b,cpx]=scmplx(A);
    besti=i;
    for j=i-1:-1:1
        tmp=A(:,[1:j-1,j+1:i,j,i+1:n+1]);
        tmp=QRZ_HH(tmp,i);
        [a,b,cpx1]=scmplx(tmp);
        if cpx1<cpx
            besti=j;
            cpx=cpx1;
        end
    end
    if besti~=i
        A=A(:,[1:besti-1,besti+1:i,besti,i+1:n+1]);
        Z=Z(:,[1:besti-1,besti+1:i,besti,i+1:n]);
        A=QRZ_HH(A,i);
    end
end
