% ALEX VERZEA
% 260324472
% COMP 540 Assignment 4, #5.
%
% createVandermondeMatrix.m    creates a Vandermonde Matrix
%
% input: m is the number of rows of A
%        n is the number of columns of A
%
% output: A is a Vandermonde Matrix
%
function A = createVandermondeMatrix(m,n)
for j = 1:n
    for i = 1:m
        A(i,j) = (j/n)^(i-1);
    end;
end;
