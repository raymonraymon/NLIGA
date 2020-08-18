function drawSphere(center, r, clr)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Draw a sphere of radius r colored clr on center
% 
% acceptable forms:
%       drawSphere(center)
%       drawSphere(center, r)
%       drawSphere(center, r, clr)
% 
% inputs:
%       center: center of the sphere, center(x,y,z)
%       r:      radius of the sphere, default = 1
%       clr:    color of the sphere, default = 'k'
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% default values
if nargin == 1
    r = 1;
    clr = 'k';
elseif nargin == 2;
    clr = 'k';
end

% center of the sphere
xc = center(1);
yc = center(2);
zc = center(3);

% division number in parallel
nTheta = 18;

% number number of meridian
nPhi = 36;

% compute spherical coordinates
theta = linspace(0, pi, nTheta);
phi = linspace(0, 2*pi, nPhi);

% convert to cartesian coordinates
x = xc + cos(phi')*sin(theta)*r;
y = yc + sin(phi')*sin(theta)*r;
z = zc + ones(length(phi),1)*cos(theta)*r;

options = {'FaceColor', clr, 'linestyle', 'none'};
surf(x, y, z, options{:});

end