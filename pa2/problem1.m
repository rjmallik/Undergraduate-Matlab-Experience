disp("Problem 1.2");
A = [0, 2, -8;
    2, -2, 1;
    -4, 5, -7];
A_aug = [A , eye(3)];
A_rref = rref(A_aug);
A_inv = A_rref(:, 4:6);
A_inv_builtin = inv(A);
% (c) Compute A*A^-1 and A^-1*A to verify that both products are I3
AA_inv = A * A_inv;
invA_A = A_inv * A;
% Display results
disp('A^-1 computed from RREF:');
disp(A_inv);
disp('A^-1 computed using built-in inv function:');
disp(A_inv_builtin);
disp('A * A^-1:');
disp(AA_inv);
disp('A^-1 * A:');
disp(invA_A);
disp(A_rref);

disp("Problem 1.3");
b2 = [0;
    0;
    0];
x_b2 = A_inv * b2;

% Display the solution set
disp('The solution set for Ax = b2 is:');
disp(x_b2);

