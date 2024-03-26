%problem1
%Part 1.1
disp('Problem 1.1');
%Part 1.1.1
disp('Part 1.1.1');
u1 = [1;
    -1;
    2;
    0];
u2 = [-3;
    1;
    -8;
    2];
u3 = [1;
    -7;
    1;
    5];
u4 = [5;
    -5;
    10;
    0];

A = [u1, u2, u3, u4];
NS = null(A);
disp('The Null Space:')
disp(NS);
disp('Null space is not empty. It is linearly dependent');

%Part 1.1.2
disp('Part 1.1.2');
B = [u1, u2, u3];
rrefB = rref(B);
disp('RREF B:');
disp(rrefB);
disp('It is Linear Independent');

%Part 1.1.3
disp('Part 1.1.3');
Dim_S = rank(B);
disp('Dimension of S:');
disp(Dim_S);

%Part 1.1.4
disp('Part 1.1.4');
v1 = [13;
    -33;
    26;
    16];

coefficient = B \ v1;
disp('Coefficient Vector Coordinates:');
disp(coefficient);

%Problem 1.2
disp('Problem 1.2');

%Part 1.2.1
disp('1.2.1');

T = [2, 0, 1, 1;
    -2, 1, -1, 0;
    0, 0, 0, -2];

rankT = rank(T);
disp('Rank at T:');
disp(rankT);

disp('Problem 1.3');

