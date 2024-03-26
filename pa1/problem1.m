disp('Problem 1.1');
B = [8, 2;
    -1, 0;
    3, 2;
    6, 5;
    9, -2];
b1 = B( : , 1);
disp('b1:');
disp(b1);
b1size = size(b1);
disp('Size of b1:');
disp(b1size);
B(:, 2) = B(:, 2) + b1;
disp('b1 added to second column of B');
disp(B)
B(3, :) = B(3, :) * 5;
disp('Multiplying third row of B by 5 and assign it to the 3rd row of B:');
disp(B);
C = B([2 3], :);
disp('Define a new variable C that contains the second and third rows of B');
disp(C);

disp('Problem 1.2');
D = [1, 2;
    -2, 4];
x1 = [1;
    0];
x2 = [0;
    1];
Dx1 = D * x1;
Dx2 = D * x2;
aggTotal = Dx1 + Dx2;
agg = D * (x1+x2);
disp('Dx1:');
disp(Dx1);
disp('Dx2:');
disp(Dx2);
disp('Dx1 + Dx2:');
disp(aggTotal);
disp('D(x1+x2):');
disp(agg);

a1 = 3;
a2 = -2;
big1 = (a1 * Dx1) + (a2 * Dx2);
disp('a1Dx1 + a2Dx2:');
disp(big1);
big2 = D*((a1*x1)+(a2*x2));
disp('D*(a1x1 + a2x2):');
disp(big2);
disp(['The relationship between the formula contianed ' ...
    'in the variable big1 and big2 is symbiotic as you get the same result ' ...
    'no matter the formula you pick']);


disp('Problem 1.3');
% Open figure #1 and make sure it is empty .
figure(1); clf
% Turn ' hold ' on so plots don ' t disappear .
hold on
%First for loop
g1 = 0;
for g2 = -5:5
    g = (g1*x1)+(g2*x2);
    plot(g1, g , 'black+');
    Dg = D*g;
    plot(g1, Dg , 'red.');
    y = (g1*Dx1) + (g2* Dx2);
    plot(g1,y , 'blueo');
end
% Second for loop
g2 = 0;
for g1 = -5:5
    g = (g1*x1)+(g2*x2);
    plot(g1, g , 'black+');
    Dg = D*g;
    plot(g1, Dg , 'red.');
    y = (g1*Dx1) + (g2* Dx2);
    plot(g1,y , 'blueo');
end
%Nested for loop
for g1 = -5:5
    for g2 = -5:5
    plot(g1, g2 , 'black+');
    end
end



