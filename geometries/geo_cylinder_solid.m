function cylinder = geo_cylinder_solid( R, L )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build cylindrical SOLID model.
% model.
% Input:
%   R - radius 
%   L - length
% 
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 2    % default parameters
    R = 1.0;
    L = 4.0;
end

% define control points and knots vector
rad = 45*pi/180;
w = cos(rad);
coefs = zeros(4,3,3,2);
coefs(:,:,1,1) = [-R*cos(rad)     0   -R*sin(rad)      1;
                           0      0   -R/cos(rad)*w    w;
                   R*cos(rad)     0    -R*sin(rad)     1]';
coefs(:,:,2,1) = [-R/cos(rad)*w   0    0      w;
                   0              0    0    1/2;
                   R/cos(rad)*w   0    0      w]';
coefs(:,:,3,1) = [-R*cos(rad)     0  R*sin(rad)     1;
                            0     0  R/cos(rad)*w   w;
                   R*cos(rad)     0  R*sin(rad)     1]';
coefs(:,:,1,2) = [-R*cos(rad)     L     -R*sin(rad)      1;
                           0      L*w   -R/cos(rad)*w    w;
                   R*cos(rad)     L     -R*sin(rad)      1]';
coefs(:,:,2,2) = [-R/cos(rad)*w   L*w    0    w;
                   0              L/2    0    1/2;
                   R/cos(rad)*w   L*w    0    w]';
coefs(:,:,3,2) = [-R*cos(rad)     L    R*sin(rad)        1;
                            0     L*w  R/cos(rad)*w      w;
                   R*cos(rad)     L    R*sin(rad)        1]';
knots{1}=[0 0 0 1 1 1];
knots{2}=[0 0 0 1 1 1];
knots{3}=[0 0 1 1];

% build nurbs solid by using control points and knots vector
cylinder = nrbmak(coefs,knots);

% degeree elevate
cylinder = nrbdegelev(cylinder,[1,2,2]);

% insert knots
RefinementX = 3;    % the number of knots inseted in u direction 
RefinementY = 3;    % the number of knots inseted in v direction
RefinementZ = 3;    % the number of knots inseted in w direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
iwknots = 1/(RefinementZ+1):1/(RefinementZ+1):RefinementZ/(RefinementZ+1);
cylinder = nrbkntins(cylinder, {iuknots,ivknots,iwknots});

% plot nurbs solid
plot_nurbs(cylinder, 0,1);
axis equal
view(3);

end

