% ALEX VERZEA
% 260324472
% COMP 540 Assignment 4, #5.
%
% classicGramSchmidt.m    Classic Gram-Schmidt Orthogonalization
%
% input: A is an m x n matrix (full column rank)
%
% output: Q is an m x n matrix
%         R is an n x n matrix
%
function [Q1, R] = classicGramSchmidt(A)
[m, n] = size(A);
Q1 = zeros(m, n);
R = zeros(n, n);
for k = 1:n,
    a = A(:,k);
    for i = 1:k-1
        R(i,k) = Q1(:,i)' * A(:,k);
    end;
    for i = 1:k-1
        a = a - R(i,k) * Q1(:,i);
    end;
    R(k,k) = norm(a);
    Q1(:,k) = a / R(k,k);
end;
