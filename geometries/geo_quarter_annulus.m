function annulus = geo_quarter_annulus(r,R)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a one-quarter annulus
% r - inner radius
% R - outer radius
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 2    % default parameters
    r = 1.0;  
    R = 2.0; 
end
% define control points and knots vector
alpha = 45/180*pi;
w = cos(alpha);
coefs = zeros(4,3,2);
coefs(:,:,1) = [0,  r, 0, 1;
                r*w,r*w,0,w;
                r,  0,  0,1]';
coefs(:,:,2) = [0,  R, 0,  1;
                R*w,R*w,0, w;
                R,  0,  0, 1]';
            
knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
annulus = nrbmak(coefs,knots);

% degeree elevate
annulus = nrbdegelev(annulus,[1 2]);

% insert knots
RefinementX = 7;    % the number of knots inseted in u direction 
RefinementY = 7;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
annulus = nrbkntins(annulus,{iuknots ivknots});

% plot nurbs patch
figure
plot_nurbs(annulus, 0, 1);
view(2);
axis equal
end

