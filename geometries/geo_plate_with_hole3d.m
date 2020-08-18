function plate_hole = geo_plate_with_hole3d
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% build a 3d square plate with hole
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% define control points and knots vector
L = 82.5;  % side length
R = 6.35;  % radius
h = 10; % height
coefs = zeros(4,4,3,2);
alpha = 22.5/180*pi;
w = cos(alpha)^2;
coefs(:,:,1,1) = [0, R, 0, 1; R*tan(alpha)*w, R*w, 0, w; R*w, R*tan(alpha)*w, 0, w; R, 0, 0, 1]';
coefs(:,:,2,1) = [0, (R+L)/2, 0, 1; (R+L)/2*tan(alpha)*w, (R+L)/2*w, 0, w; (R+L)/2*w, (R+L)/2*tan(alpha)*w, 0, w; (R+L)/2, 0, 0, 1]';
coefs(:,:,3,1) = [0, L, 0, 1; L*w, L*w, 0, w; L*w, L*w, 0, w; L, 0, 0, 1]';
coefs(:,:,1,2) = [0, R, h, 1; R*tan(alpha)*w, R*w, h*w, w; R*w, R*tan(alpha)*w, h*w, w; R, 0, h, 1]';
coefs(:,:,2,2) = [0, (R+L)/2, h, 1; (R+L)/2*tan(alpha)*w, (R+L)/2*w, h*w, w; (R+L)/2*w, (R+L)/2*tan(alpha)*w, h*w, w; (R+L)/2, 0, h, 1]';
coefs(:,:,3,2) = [0, L, h, 1; L*w, L*w, h*w, w; L*w, L*w, h*w, w; L, 0, h, 1]';
knots{1} = [0 0 0 0.5 1 1 1];
knots{2} = [0 0 0 1 1 1];
knots{3} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
plate_hole = nrbmak(coefs, knots);

% degeree elevate
plate_hole = nrbdegelev(plate_hole,[1,1,2]);

% insert knots
RefinementX = 4;    % the number of knots inseted in u direction 
RefinementY = 3;    % the number of knots inseted in v direction
RefinementZ = 0;    % the number of knots inseted in w direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
iwknots = 1/(RefinementZ+1):1/(RefinementZ+1):RefinementZ/(RefinementZ+1);
plate_hole = nrbkntins(plate_hole, {iuknots ivknots iwknots});

% plot nurbs solid
plot_nurbs(plate_hole, 0,1);
axis equal
end

