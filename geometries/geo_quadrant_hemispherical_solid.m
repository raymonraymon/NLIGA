function hemisphere = geo_quadrant_hemispherical_solid( R, phi, t )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a quadrant of hemispherical SOLID model.
% Note that 'geo_quadrant_hemisphere_shell' is the middle surface of this
% model when given the same input parameters.
% model.
% Input:
%   R - radius 
%   phi - angle
%   t - thickness
% 
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 3    % default parameters
    R = 10;
    phi = 18;
    t = 0.04;
end

% define control points and knots vector
rad1 = phi/180*pi;
rad2 = (90-phi)/2/180*pi;
w1   = cos(45/180*pi);
w2   = cos((90-phi)/2/180*pi);
r1   = R - t/2;
r2   = R + t/2;
coefs = zeros(4,3,3,2);
coefs(:,:,1,1) = [r1, 0, 0, 1; r1*w1, r1*w1, 0, w1; 0, r1, 0, 1; ]';
coefs(:,:,2,1) = [   r1*w2,       0,     r1*tan(rad2)*w2,    w2; 
                  r1*w1*w2, r1*w1*w2, r1*tan(rad2)*w1*w2, w1*w2;
                        0,     r1*w2,    r1*tan(rad2)*w2,    w2;]';
coefs(:,:,3,1) = [   r1*sin(rad1),               0,  r1*cos(rad1),     1; 
                  r1*sin(rad1)*w1, r1*sin(rad1)*w1,  r1*cos(rad1)*w1, w1; 
                                0,    r1*sin(rad1),  r1*cos(rad1),     1;]';
coefs(:,:,1,2) = [r2, 0, 0, 1; r2*w1, r2*w1, 0, w1; 0, r2, 0, 1; ]';
coefs(:,:,2,2) = [   r2*w2,       0,     r2*tan(rad2)*w2,    w2; 
                  r2*w1*w2, r2*w1*w2, r2*tan(rad2)*w1*w2, w1*w2;
                        0,     r2*w2,    r2*tan(rad2)*w2,    w2;]';
coefs(:,:,3,2) = [   r2*sin(rad1),               0,  r2*cos(rad1),     1; 
                  r2*sin(rad1)*w1, r2*sin(rad1)*w1,  r2*cos(rad1)*w1, w1; 
                                0,    r2*sin(rad1),  r2*cos(rad1),     1;]';
knots{1}=[0 0 0 1 1 1];
knots{2}=[0 0 0 1 1 1];
knots{3}=[0 0 1 1];

% build nurbs solid by using control points and knots vector
hemisphere = nrbmak(coefs,knots);

% degeree elevate
hemisphere = nrbdegelev(hemisphere,[1,1,2]);

% insert knots
RefinementX = 7;    % the number of knots inseted in u direction 
RefinementY = RefinementX;    % the number of knots inseted in v direction
RefinementZ = 0;    % the number of knots inseted in w direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
iwknots = 1/(RefinementZ+1):1/(RefinementZ+1):RefinementZ/(RefinementZ+1);
hemisphere = nrbkntins(hemisphere, {iuknots,ivknots,iwknots});

% plot nurbs solid
plot_nurbs(hemisphere, 0,1);
axis equal
view(3);

end

