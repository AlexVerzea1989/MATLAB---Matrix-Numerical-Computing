fprintf('1. For m = 20 and n = 10');
A = createVandermondeMatrix(20,10);
fprintf('\n');

fprintf('Classic Gram Schmidt');
fprintf('\n');
[Q1, R] = classicGramSchmidt(A);
fprintf('\n');
fprintf('norm(A-Q1*R)/norm(A)');
norm(A-Q1*R)/norm(A)
fprintf('\n');
fprintf('norm((Q1^T)*Q1-I)');
norm((Q1'*Q1)-eye(10))
norm(Q1')
fprintf('\n');

fprintf('Modified Gram Schmidt');
fprintf('\n');
[Q1, R] = modifiedGramSchmidt(A);
fprintf('\n');
fprintf('norm(A-Q1*R)/norm(A)');
norm(A-Q1*R)/norm(A)
fprintf('\n');
fprintf('norm((Q1^T)*Q1-I)');
norm(Q1'*Q1-eye(10))
fprintf('\n');

fprintf('qr(A)');
fprintf('\n');
[Q1,R] = qr(A);
fprintf('\n');
fprintf('Note that QR Factorization is implemented a bit differently in MATLAB. [Q,R] = qr(A), where A is m-by-n, produces an m-by-n upper triangular matrix R and an m-by-m unitary matrix Q so that A = Q*R.');
fprintf('\n');
fprintf('\n');
fprintf('norm(A-Q1*R)/norm(A)');
norm(A-Q1*R)/norm(A)
fprintf('\n');
fprintf('norm((Q1^T)*Q1-I)');
norm(Q1'*Q1-eye(20))
fprintf('\n');

fprintf('2. For m = 30 and n = 20');
A = createVandermondeMatrix(30,20);
fprintf('\n');

fprintf('Classic Gram Schmidt');
fprintf('\n');
[Q1, R] = classicGramSchmidt(A);
fprintf('\n');
fprintf('norm(A-Q1*R)/norm(A)');
norm(A-Q1*R)/norm(A)
fprintf('\n');
fprintf('norm((Q1^T)*Q1-I)');
norm(Q1'*Q1-eye(20))
fprintf('\n');

fprintf('Modified Gram Schmidt');
fprintf('\n');
[Q1, R] = modifiedGramSchmidt(A);
fprintf('\n');
fprintf('norm(A-Q1*R)/norm(A)');
norm(A-Q1*R)/norm(A)
fprintf('\n');
fprintf('norm((Q1^T)*Q1-I)');
norm(Q1'*Q1-eye(20))
fprintf('\n');

fprintf('qr(A)');
fprintf('\n');
[Q1,R] = qr(A);
fprintf('\n');
fprintf('Note that QR Factorization is implemented a bit differently in MATLAB. [Q,R] = qr(A), where A is m-by-n, produces an m-by-n upper triangular matrix R and an m-by-m unitary matrix Q so that A = Q*R.');
fprintf('\n');
fprintf('\n');
fprintf('norm(A-Q1*R)/norm(A)');
norm(A-Q1*R)/norm(A)
fprintf('\n');
fprintf('norm((Q1^T)*Q1-I)');
norm(Q1'*Q1-eye(30))
fprintf('\n');
