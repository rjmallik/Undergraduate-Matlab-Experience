disp("Problem 2.1");
% Define the matrix A
A = [0 3 3; 
    -4 5 2; 
    8 -4 3];

% Perform PLU factorization
% Define a sample matrix A
% Get the size of A
[n, ~] = size(A);

% Initialize P, L, and U
P = eye(n);
L = eye(n);
U = A;

for k = 1:n-1
    % Partial Pivoting
    [~, maxIndex] = max(abs(U(k:n, k)));
    maxIndex = maxIndex + k - 1;

    % Swap rows in U
    temp = U(k, :);
    U(k, :) = U(maxIndex, :);
    U(maxIndex, :) = temp;

    % Swap rows in P
    temp = P(k, :);
    P(k, :) = P(maxIndex, :);
    P(maxIndex, :) = temp;

    % Swap rows in L
    if k > 1
        temp = L(k, 1:k-1);
        L(k, 1:k-1) = L(maxIndex, 1:k-1);
        L(maxIndex, 1:k-1) = temp;
    end

    % Elimination
    for j = k+1:n
        L(j, k) = U(j, k) / U(k, k);
        U(j, :) = U(j, :) - L(j, k) * U(k, :);
    end
end

% Compute PA and LU
PA = P * A;
LU = L * U;

% Display the results
disp('Matrix P:');
disp(P);
disp('Matrix L:');
disp(L);
disp('Matrix U:');
disp(U);
disp('Product PA:');
disp(PA);
disp('Product LU:');
disp(LU);


% Define the vector b
b = [-18; -32; 22];

% Solve Ax = b using the PLU factorization
y = L\(P * b);
x = U\y;

% Verify that Ax = b
Ax = A*x;

% Display the solution x and the product Ax
disp('The solution x:');
disp(x);
disp('The product Ax:');
disp(Ax);

disp("Problem 2.3");
[L ,U , P ] = lu ( A );
% Given matrix A from the previous problem
A = [0, 3, 3; -4, 5, 2; 8, -4, 3];

% Given vector b, which you will need to define
b = [-18;
    -32;
    22]; % Replace with the actual values

% Perform LU decomposition using MATLAB's built-in function
[L, U, P] = lu(A);

% Solve for y in Ly = Pb
y = L \ (P*b);

% Display the result for y
disp('The solution vector y is:');
disp(y);

% Solve for x in Ux = y
x = U \ y;

% Display the result for x
disp('The solution vector x is:');
disp(x);
