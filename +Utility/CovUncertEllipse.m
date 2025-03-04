function fig = CovUncertEllipse(finalposCov, scRefpos, xhist, fig)
% Eigen decomposition: V columns are eigenvectors
[eigVec, eigVal] = eig(finalposCov);

% The mean position (estimated state + final reference position)
Xfinal = xhist(1:3, end) + scRefpos(end, :)';  % Assuming Xfinal is 3x1 vector

% Get the eigenvalues (axis lengths of the ellipse) and eigenvectors
axis_length = 6*sqrt(diag(eigVal));  % Length of the axes of the ellipse
theta = linspace(0, 2*pi, 100);  % Angles for ellipse points

% Parametric equations for the ellipse
ellipse_points = [cos(theta); sin(theta)];
ellipse_points = [ellipse_points(1,:) * axis_length(1);
    ellipse_points(2,:) * axis_length(2)];

% get each sigma value
sigma_x = num2str(axis_length(1), '%1.d');
sigma_y = num2str(axis_length(2), '%1.d');
sigma_z = num2str(axis_length(3), '%1.d');

% -- XY Covariance Plot
figure(fig);
% Rotate the ellipse using the first two eigenvectors (for XY plane)
rotated_ellipse_xy = eigVec(:, 1:2) * ellipse_points;
plot(Xfinal(1) + rotated_ellipse_xy(1,:), Xfinal(2) + rotated_ellipse_xy(2,:), 'r', 'LineWidth', 2);
hold on
plot(Xfinal(1), Xfinal(2), '.', 'LineWidth', 50)
xlabel('X [km]');
ylabel('Y [km]');
title('XY Plane Covariance Ellipse');
grid on;
hold off;

% Create the annotation string
str = {['\sigma_x = ' sigma_x], ...
    ['\sigma_y = ' sigma_y]};

% Create annotation box in the top-right corner
annotation('textbox', [0.7, 0.75, 0.15, 0.15], 'String', str, 'BackgroundColor', 'white');

fig = fig + 1;

% -- XZ Covariance Plot
figure(fig);
% Rotate the ellipse using the first and third eigenvectors (for XZ plane)
rotated_ellipse_xz = eigVec(:, [1, 3]) * ellipse_points;
plot(Xfinal(1) + rotated_ellipse_xz(1,:), Xfinal(3) + rotated_ellipse_xz(2,:), 'r', 'LineWidth', 2);
hold on
plot(Xfinal(1), Xfinal(3), '.', 'LineWidth', 50)
axis equal;
xlabel('X [km]');
ylabel('Z [km]');
title('XZ Plane Covariance Ellipse');
grid on;
hold off;

% Create the annotation string
str = {['\sigma_x = ' sigma_x], ...
    ['\sigma_z = ' sigma_z]};

% Create annotation box in the top-right corner
annotation('textbox', [0.7, 0.75, 0.15, 0.15], 'String', str, 'BackgroundColor', 'white');


fig = fig + 1;

% -- YZ Covariance Plot
figure(fig);
% Rotate the ellipse using the second and third eigenvectors (for YZ plane)
rotated_ellipse_yz = eigVec(:, 2:3) * ellipse_points;
plot(Xfinal(2) + rotated_ellipse_yz(1,:), Xfinal(3) + rotated_ellipse_yz(2,:), 'r', 'LineWidth', 2);
hold on
plot(Xfinal(2), Xfinal(3),'.', 'LineWidth', 50)
axis equal;
xlabel('Y [km]');
ylabel('Z [km]');
title('YZ Plane Covariance Ellipse');
grid on;
hold off;

% Create the annotation string
str = {['\sigma_y = ' sigma_y], ...
    ['\sigma_z = ' sigma_z]};

% Create annotation box in the top-right corner
annotation('textbox', [0.7, 0.75, 0.15, 0.15], 'String', str, 'BackgroundColor', 'white');

end
