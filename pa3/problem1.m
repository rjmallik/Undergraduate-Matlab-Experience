disp("Problem 1.1:");
%1.1.1
A = [1 -3 1 5; -1 1 -7 -5; 2 -8 1 10; 0 2 5 0];
rank_A = rank(A);
disp(['Rank of A: ', num2str(rank_A)]);
isDependent = (rank_A < size(A, 2));
disp(['Is the set linearly dependent(1 is true anything else is false): ', num2str(isDependent)]);
null_space_A = null(A, 'r');
c = null_space_A(:, end); % Coefficients c1, c2, c3
disp(['Coefficients c1, c2, c3: ', num2str(c')]);
%1.1.2
B = [1 -3 1; -1 1 -7; 2 -8 1; 0 2 5];
rank_B = rank(B);
disp(['Rank of B: ', num2str(rank_B)]);
isIndependent = (rank_B == size(B, 2));
disp(['Are the vectors linearly independent(1 is true anything else is false): ', num2str(isIndependent)]);
%1.1.4
% Given vectors u1, u2, u3, and v
u1 = [1; -1; 2; 0];
u2 = [-3; 1; -8; 2];
u3 = [1; -7; 1; 5];
u4 = [5; -5; 10; 0]; % Redundant vector
v = [13; -33; 26; 16];

% Form matrix B with u1, u2, u3
B = [u1 u2 u3];

% Solve for the coefficients c in Bc = v
c = B \ v;

% Display the coordinate vector [v]_B
disp('Coordinate Vector of [v] relative to B:');
disp(c');
% Calculate the linear combination using the obtained coefficients(To check
% answers and see if it is outputs the same as v.)
linear_combination = B * c;

disp(['Linear Combination using coefficients: ', num2str(linear_combination')]);
% Check if the linear combination equals the original vector v

disp("Problem 1.2:");
%1.2.4 & 1.2.5
A_T = [2 0 1 1; -2 1 -1 0; 0 0 0 -2];
Tref = rref(A_T);
disp(Tref); 
% Calculate the null space of A_T
null_space_T = null(A_T, 'r');

% Find the dimension of the null space
nullity_T = size(null_space_T, 2);

disp(['Nullity of A_T(If the nullity_T = 0, then T is one-to-one): ', num2str(nullity_T)]);