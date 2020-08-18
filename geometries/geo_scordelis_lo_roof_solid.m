function roof = geo_scordelis_lo_roof_solid( L, R, t, theta )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build Scordelis-Lo roof SOLID model.
% Note that 'geo_scordelis_lo_roof_shell' is the middle surface of this
% model.
% Input:
%   R - radius 
%   L - length
%   t - thickness
%   theta - angle
% 
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 4    % default parameters
    R = 25;
    L = 50;
    t = 0.25;
    theta = 40;
end
rad = theta/180*pi;
w = cos(rad);
r1 = R - t/2;
r2 = R + t/2;
% define control points and knots vector
coefs = zeros(4,3,2,2);
coefs(:,:,1,1) = [-r1*sin(rad), 0, r1*cos(rad), 1;
                0, 0, r1/cos(rad)*w, w;
                r1*sin(rad), 0, r1*cos(rad), 1;]';
coefs(:,:,2,1) = [-r1*sin(rad), L, r1*cos(rad), 1;
                0, L*w, r1/cos(rad)*w, w;
                r1*sin(rad), L, r1*cos(rad), 1;]';
coefs(:,:,1,2) = [-r2*sin(rad), 0, r2*cos(rad), 1;
                0, 0, r2/cos(rad)*w, w;
                r2*sin(rad), 0, r2*cos(rad), 1;]';
coefs(:,:,2,2) = [-r2*sin(rad), L, r2*cos(rad), 1;
                0, L*w, r2/cos(rad)*w, w;
                r2*sin(rad), L, r2*cos(rad), 1;]';
knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 1 1];
knots{3} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
roof = nrbmak(coefs,knots);

% degeree elevate
roof = nrbdegelev(roof,[1,2,1]);

% insert knots
RefinementX = 7;    % the number of knots inseted in u direction 
RefinementY = RefinementX;    % the number of knots inseted in v direction
RefinementZ = 0;    % the number of knots inseted in w direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
iwknots = 1/(RefinementZ+1):1/(RefinementZ+1):RefinementZ/(RefinementZ+1);
roof = nrbkntins(roof, {iuknots,ivknots,iwknots});

% plot nurbs solid

plot_nurbs(roof, 0,1);

view(3);
axis equal


end

