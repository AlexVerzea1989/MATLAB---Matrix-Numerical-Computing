% ALEX VERZEA
% 260324472
% COMP 540 A6 4
% 
% Factorization, recombination, iteration. No stopping criteria used.
%
% For matrices that converge quickly, I put a prompt line so that you 
% can stop between each iteration to check the validity of the results.
% For matrices converging slowly, I used a while iteration < (big number)
% mechanism to show the matrix at each iteration.
%
% eig(A) will print the eigenvalues both for the starting matrix and for
% each iteration. I have put it in comments here because for the matrices
% requiring many iterations, the eigenvalues take a lot of space on the
% screen.
%
% Input: Use one of these exact inputs.
% A = [[5 4 1 1]' [4 5 1 1]' [1 1 4 2]' [1 1 2 4]'];
% A = [[6 4 4 1]' [4 6 1 4]' [4 1 6 4]' [1 4 4 6]'];
% A = [[33 -24 -8]' [16 -10 -4]' [72 -57 -17]'];
% A = [[6 4 4 4]' [-3 2 -2 2]' [4 4 3 3]' [1 0 1 1]'];
% A = [[4 0 5 3]' [-5 4 -3 0]' [0 -3 4 5]' [3 -5 0 4]'];
% A = [[10 9 8 6 4 2]' [-19 -18 -16 -12 -8 -4]' [17 17 15 12 8 4]'
% [-12 -12 -11 -10 -6 -3]' [4 4 4 4 1 1]' [1 1 1 1 2 0]']
%
% Return: A
%
function I = recombination(A)
%eig(A)

iteration = 1;
while iteration < 30001
    iteration
    [Q R] = qr(A);
    A = R * Q
    %eig(A)
    iteration = iteration + 1;
    %input('Press Enter for next iteration.')
end

end
