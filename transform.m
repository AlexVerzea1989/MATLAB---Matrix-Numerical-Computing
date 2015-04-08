% ALEX VERZEA
% 260324472
% COMP 540 A5 4a
% 
% First we need to make sure A is symmetric. Then we find a Jacobi
% Rotation to make the matrix diagonal.
%
% Input: Matrix A 2x2
% Return: Matrices J1 and J2 which multiplied with A give the diagonal
% matrix on the right([d1 0] [0 d2])
%
function [J1,J2]=transform(A)
B=[1 0;0 1];

% If A is not symmetric.
if A(1,2) ~= A(2,1) 
    if (A(1,1) + A(2,2)) == 0
        c = 0;
        s = 1;
    else
        t = (A(1,2) - A(2,1)) / (A(1,1) + A(2,2));
        c = 1 / sqrt(1 + t^2);
        s = t * c;
    end
    B = [c -s ; s c];
end
A = B * A;

% Now, we find a Jacobi Rotation.
if A(1,2) == 0
    c1 = 1;
    s1 = 0;
else
    p = (A(2,2) - A(1,1)) / (2 * A(1,2));
    t1 = -p -sign(p) * sqrt(p^2 + 1);
    t2 = -1 / t1;
    if abs(t1) > abs(t2)
        t = t2;
    else
        t = t1;
    end
    c1 = 1 / sqrt(1 + t^2);
    s1 = t * c1;
end

% Put it all together. J1 and J2 are the two matrices we need for c's and
% s's.
B1 = [c1 -s1 ; s1 c1];
B2 = [c1 s1 ; -s1 c1];
J1 = B * B1;
J2 = B2;
end
