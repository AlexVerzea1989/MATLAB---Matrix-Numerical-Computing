function srate=srate_babai(R,sigma,n)
srate = 1;
if nargin<3
    n=size(R,2);
    m=size(R,1);
    if m<=n
        n=m;
    else
        n=n-1;
    end
end
for i=1:n
    srate=srate*(2*normcdf(abs(R(i,i))/(2*sigma))-1);
end
