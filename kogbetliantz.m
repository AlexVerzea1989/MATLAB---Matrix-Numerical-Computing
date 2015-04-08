% ALEX VERZEA
% 260324472
% COMP 540 A5 4b
% 
% Using 4a repeatedly to find the SVD of a given real m x n matrix.
% m, n >= 0. This is Kogbeliantz from the course notes.
%
% Input: Matrix A 2x2
% Return: SVD decomposition of A.
%
function [U,A,V]=kogbetliantz(A)
[m,n] = size(A);
U = eye(m);
V = eye(n);
e = 0.000000001 * norm(A,'fro');

while (sum(sum(A.^2)) - sum(diag(A).^2)) > e
    for i = 1 : n
        for j = i + 1 : n
            Ap = [A(i,i) A(i,j) ; A(j,i) A(j,j)];
            [J1,J2] = transform(Ap);
            J1p = eye(m);
            J2p = eye(n);
            J1p(i,i) = J1(1,1);
            J1p(i,j) = J1(1,2);
            J1p(j,i) = J1(2,1);
            J1p(j,j) = J1(2,2);
            J2p(i,i) = J2(1,1);
            J2p(i,j) = J2(1,2);
            J2p(j,i) = J2(2,1);
            J2p(j,j) = J2(2,2);
            A = J1p * A * J2p;
            U = U * transpose(J1p);
            V = V * J2p;
        end
        for j = n+1 : m
            Ap = [A(i,i) A(j,i) ; A(j,i) A(j,i)];
            [J1,J2] = transform(Ap);
            J1p = eye(m);
            J1p(i,i) = J1(1,1);
            J1p(i,j) = J1(1,2);
            J1p(j,i) = J1(2,1);
            J1p(j,j) = J1(2,2);
            A = J1p * A;
            U = U * transpose(J1p);
        end
    end
end

% Singular values must be > 0
E = eye(m);
for i = 1 : n
    if A(i,i) < 0
        E(i,i) = -1;
    end
end
A = E * A;

% U = E * U;
% V =E(1:n , 1:n) * V;
% Must have a11 > a22 > ...
% And if we don't have it, we pivot until we do.
for i = 1 : n-1
    for j = i+1 : n
        if A(i,i) < A(j,j)
            A([i,j],:) = A([j,i],:);
            A(:,[i,j]) = A(:,[j,i]);
            U(:,[i,j]) = U(:,[j,i]);
            V(:,[i,j]) = V(:,[j,i]);
        end
    end
end
end
