
%%%%%%%%%%%%%%%%%
% Problem Setup %
%%%%%%%%%%%%%%%%%

clear all
load('problem2_data.mat')

% The height of the image is the number of rows in each array.
img_height = size(X_pos, 1); % Equals size(Y_pos, 1)

% The width of the image is the number of rows in each array.
img_width = size(X_pos, 2); % Equals size(Y_pos, 2) 

% Open and clear figure 1. 
figure(1); clf

% Plot the image.
Z_pos = zeros(img_height, img_width);
surf(X_pos, Y_pos, Z_pos, RGB, 'edgecolor', 'none')

axis tight
axis equal
box on
grid off
view(0, 90)

title('What is this?')

%%
%%%%%%%%%%%%%%%
% Problem 2.1 %
%%%%%%%%%%%%%%%

% TODO: Define A1, A2, A3
A1 = [cos(5*pi/6), -sin(5*pi/6); sin(5*pi/6), cos(5*pi/6)];
A2 = [1, 1; 0, 1];
A3 = [2, 0; 0, 1/4];

c_ll = [-1; -1]; % Lower-left corner
c_lr = [ 1; -1]; % Lower-right corner
c_ur = [ 1;  1]; % Upper-right corner
c_ul = [-1;  1]; % Upper-left corner

% Create a 2x5 array containing the position of each corner in its columns.
Sq = [c_ll, c_lr, c_ur, c_ul, c_ll];

% Apply a transformation to the square (moving each corner) and store it in the
% next index of the array.
Sq(:,:,2) = A1*Sq(:,:,1);
Sq(:,:,3) = A2*Sq(:,:,2);
Sq(:, :, 4) = A3 * Sq(:, :, 3);

fig2 = figure(2); clf
hold on
box on
axis equal

for k = 1:size(Sq, 3)
    plot(Sq(1,:,k), Sq(2,:,k), '*-', 'linewidth', 2);
    hold on; % Add this line to overlay plots
end

title('Evolution of a square')
xlabel('x')
ylabel('y')
legend('Original', 'Rotate (A_1)', ...
'Shear (A_2 A_1)', 'Scale (A_3 A_2 A_1)');

ax2 = gca;
ax2.XAxisLocation = 'origin';
ax2.YAxisLocation = 'origin';

% Save the figure. You can change the file extension to
% other formats such as '.png' if a different format is easier
% to include in your report.
saveas(gcf, 'squares.pdf')

%%
%%%%%%%%%%%%%%%
% Problem 2.2 % 
%%%%%%%%%%%%%%%

XY_data = [reshape(X_pos, 1, numel(X_pos));
            reshape(Y_pos, 1, numel(Y_pos))];

% Solves for original coordinates

% Assuming A is the matrix representing the composition of transformations
A = A3 * A2 * A1;

% Calculate the inverse of A
A_inv = inv(A);

% XY_data is a matrix where each column represents a vector x_data
XY_original = A \ XY_data;

figure (3) ; clf
% Convert back from a two row array to arrays in the shape of the image. 
X_original = reshape(XY_original(1,:), img_height, img_width);
Y_original = reshape(XY_original(2, :), img_height, img_width);

Z_pos = zeros ( img_height , img_width ) ;
surf(X_original, Y_original, Z_pos, RGB, 'EdgeColor', 'none', 'FaceColor', 'interp');
% Prettify the plot .
axis tight
axis equal
box on
grid off
view (0 , 90)
title('Corrected Picture');
% Save the figure . You can change the file extension to
% other formats such as '. png ' if a different format is
% to include in your report .
saveas ( gcf , ' corrected_picture.png ')

Displaying problem3_template.m.