% ALEX VERZEA
% 260324472
% COMP 540 A5 4e
% 
% Moore-Penrose Pseudoinverse is A+ = V E+ U^T
%
% Input: Matrix A 2x2
% Return: A+
%
function I = moorepenrose(A)
[U,S,V] = kogbetliantz(A);
[m,n] = size(S);

% Need index of non-zero singular value.
c = 1;
for i = 1:n
    if S(i,i) > 0.0001, c = i; end
end
S(1 : c , 1 : c) = inv(S(1 : c, 1 : c));
I = V * transpose(S) * transpose(U);
end
