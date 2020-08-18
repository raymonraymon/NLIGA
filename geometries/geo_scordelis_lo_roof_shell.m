function roof = geo_scordelis_lo_roof_shell( L, R, theta )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build Scordelis-Lo roof SHELL model
% Input:
%   R - radius 
%   L - length
%   theta - angle
% 
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 3    % default parameters
    R = 25;
    L = 50;
    theta = 40;
end
rad = theta/180*pi;
w = cos(rad);

% define control points and knots vector
coefs = zeros(4,3,2);
coefs(:,:,1) = [-R*sin(rad), 0, R*cos(rad), 1;
                0, 0, R/cos(rad)*w, w;
                R*sin(rad), 0, R*cos(rad), 1;]';
coefs(:,:,2) = [-R*sin(rad), L, R*cos(rad), 1;
                0, L*w, R/cos(rad)*w, w;
                R*sin(rad), L, R*cos(rad), 1;]';
knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
roof = nrbmak(coefs,knots);

% degeree elevate
roof = nrbdegelev(roof,[1,2]);

% insert knots
RefinementX = 11;    % the number of knots inseted in u direction 
RefinementY = RefinementX;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
roof = nrbkntins(roof, {iuknots,ivknots});

% plot nurbs solid

plot_nurbs(roof, 0,1);

view(3);
axis equal


end

