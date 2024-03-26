disp('Problem 2.3');
A = [300, 900, 0;
    1, 1, 1;
    500, 2000, 0];
b = [30000;
    60;
    65000];
Ab = [A, b];
Ab_original = Ab;
Ab([1, 3], :) = Ab([3, 1], :);
Ab (1 ,:) = Ab (1 ,:) / 500;
Ab (2 ,:) = Ab(2, :) - Ab(1, :);
store1 =Ab(1, :)*300;
Ab (3, :) = Ab(3, :) - store1;
Ab([2, 3], :) = Ab([3, 2], :);
Ab (2 ,:) = Ab (2 ,:) / -300;
store2 =Ab(2, :)*4;
Ab (1, :) = Ab(1, :) - store2;
store3 = Ab(2, :) * (-3);
Ab(3, :) = Ab(3, :) - store3;
disp('Ab in its RREF:')
disp(Ab);

disp('Problem 2.4');
Ab_rref = rref(Ab_original);

% To make sure code above which calculates RREF of that specific matrix is
% right
%if isequal(Ab_rref, Ab)
%    disp('The code to find the RREF of AB is right');
%else
%    disp('The RREF of Ab is not right.');
%end

disp('RREF of Ab:');
disp(Ab_rref);
AbRank = rank(Ab);
disp('Rank of Ab:')
disp(AbRank);
colmnNum = size(Ab ,2);
disp('Number of columns:')
disp(colmnNum);
disp('Nullity of Ab:')
nullity = colmnNum - AbRank;
disp(nullity);
% nullity = columns - rank(number of pivots in the 1's though RREF)

disp('Problem 2.5');
x_sol = A \ b;
disp(x_sol);
x_star = [10; 
    30;
    20]; 
% Replace with the actual values you calculated by hand
if isequal(x_sol, x_star)
    disp('The solution matches the one calculated by hand.');
else
    disp('The solution does not match. Check your calculations.');
end

disp('Problem 2.6 - Extra Credit');
newA = [300, 900, 0, 2100;
    1, 1, 1, 1;
    500, 2000, 0, 5000];
newAb = [newA, b];
refofnewAB = rref(newAb);
disp(refofnewAB);