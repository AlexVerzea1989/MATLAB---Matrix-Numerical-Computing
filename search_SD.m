function [Zhat,rsd] = search_SD(R,y,r)
%
% Zhat = search(R,y,p) produces p optimal solutions to the upper triangular
%        integer least squares problem min_{z}||y-Rz|| by a search algorithm.
%
% Input arguments:
%    R ---- n by n real nonsingular upper triangular matrix
%    y ---- n-dimensional real vector
%    p ---- the number of optimal solutions and its default value is 1
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

while 1
        %-----------
        %fprintf('1');
        if k~=1
            fprintf('\t%c',ones(k-1,1)*'-');
        end
        fprintf('\t%d',z(k:n));
        fprintf('\n');
        %----------
    newprsd = prsd(k) + gamma*gamma;
    if newprsd < beta
        if k ~= 1 % move to level k-1
            S(1:k,k) = R(1:k,k)*z(k) + S(1:k,k+1);
            k = k - 1;
            prsd(k) = newprsd;
            c(k) = (y(k)-S(k,k+1))/R(k,k);
            z(k) = round(c(k));
            gamma = R(k,k)*(c(k)-z(k));
            if c(k) > z(k)
                d(k) = 1;
            else
                d(k) = -1;
            end
        else % a new point is found, update the set of candidate solutions
            if r>=0
                Zhat = z;
                beta = newprsd;
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
        if k == n   % the p optimal solutions have been found
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
if r>=0
    rsd = beta;
end
