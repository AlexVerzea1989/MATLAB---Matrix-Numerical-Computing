function [Zhat,rsd] = search_NSD(R,y,r)
%
% Zhat = search(R,y,p) produces p optimal solutions to the upper triangular
%        integer least squares problem min_{z}||y-Rz|| by a search algorithm.
% 
% Input arguments:
%    R ---- n by n real nonsingular upper triangular matrix
%    y ---- n-dimensional real vector
%    r ---- r>0 shrink. r<0 no shrink. r=-inf, shrink at the first hit.
%
% Output arguments:
%    Zhat - n by p integer matrix (in double precision). Its j-th column
%           s the j-th optimal solution, i.e., its residual norm is the j-th
%           smallest, so ||y-R*Zhat(:,1)|| <= ...<= ||y-R*Zhat(:,p)||


% Check input arguments
if nargin < 2 % input error
    error('Not enough input arguments!')
end

if nargin < 3
    r = inf;
end

[m,n] = size(R);
[n2,n3] = size(y);
if m ~= n || n ~= n2 || n3 ~= 1  % input error
    error('Input arguments have a matrix dimension error!')
end

% Initialization
z = zeros(n,1); % the current point
c = zeros(n,1); % c(k)=(y(k)-R(k,k+1:n)z(k+1:n))/R(k,k)
d = zeros(n,1); % d(k) is used to compute z(k)
prsd = zeros(n,1); % partial squared residual norm for z
% prsd(k)=norm(y(k+1:n)-R(k+1:n,k+1:n)z(k+1:n))^2
S = zeros(n,n+1); % S(:,k) = R(:,k:n)*z(k:n), k=1:n
Zhat = zeros(n,1); % the p candidate solutions (or points)
rsd = inf; % squared residual norms of the p candidate solutions

beta = abs(r);  % the initial ellipsoid bound


c(n) = y(n)/R(n,n);
z(n) = round(c(n));
gamma = R(n,n)*(c(n)-z(n));
if c(n) > z(n)
    d(n) = 1;
else
    d(n) = -1;
end

k = n;
i=1;
while 1
    if i>100000000
        disp(Zhat)
        keyboard;
    end
    i=i+1;
    newprsd = prsd(k) + gamma*gamma;
    if newprsd < beta
        if k ~= 1 % move to level k-1
            k1=k-1; %1
            sk=R(k1,k)*z(k)+S(k1,k); %2 extra
            c(k1)=(y(k1)-sk)/R(k1,k1); %2
            z(k1)=round(c(k1)); %1
            prsd(k1)=newprsd;
            gammak=R(k1,k1)*(c(k1)-z(k1)); %2
            if(gammak^2 +newprsd < beta) %need to go down
                %go down
                gamma=gammak;
                k=k1;
                S(1:k,k) = R(1:k,k+1)*z(k+1) + S(1:k,k+1);% TODO not optimal % O(2k)
                if c(k) > z(k)
                    d(k) = 1;
                else
                    d(k) = -1;
                end
            else % don't actually go down
                gamma=gammak;
                k=k1;
            end
        else % a new point is found, update the set of candidate solutions
            if r>=0 || isinf(beta)
                Zhat = z;
                beta = newprsd;
                rsd=beta;
                k=2;
            elseif newprsd < rsd
                rsd = newprsd;
                Zhat = z;
            end
            z(k) = z(k) + d(k);
            gamma = R(k,k)*(c(k)-z(k));
            if d(k) > 0
                d(k) = -d(k) - 1;
            else
                d(k) = -d(k) + 1;
            end
        end
    else
        if k == n   % the optimal solutions have been found
            break
        else  % move back to level k+1
            k = k + 1;
            z(k) = z(k) + d(k);
            gamma = R(k,k)*(c(k)-z(k));
            if d(k) > 0
                d(k) = -d(k) - 1;
            else
                d(k) = -d(k) + 1;
            end
        end
    end
end
