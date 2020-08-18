function drawCylinder(p1, p2, r, clr)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Draw a cylinder of radius r colored clr between p1 and p2
% 
% acceptable forms:
%       drawCylinder(p1, p2)
%       drawCylinder(p1, p2, r)
%       drawCylinder(p1, p2, r, clr)
% 
% inputs:
%       p1,p2: centers of the bottom of the cylinder, p1(x1,y1,z1), p2(x2,y2,z2)
%       r:     radius of the cylinder, default = 1
%       clr:   color of the cylinder, default = 'w'
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% default values
if nargin == 2
    r = 1;
    clr = 'w';
elseif nargin == 3;
    clr = 'w';
end

% division number in circumference
nPhi = 36;
phi = linspace(0, 2*pi, nPhi);

% relative coordinate system in the bottom of the cylinder
n = p1-p2;
n = n/sqrt(sum(n.^2));
a = zeros(1,3);
a(1) = p2(1) + n(1)*cos(pi/4);
a(2) = p2(2) + n(2)*sin(pi/4);
a(3) = 1;
r1 = cross(n,a);
r1 = r1/sqrt(sum(r1.^2))*r;
r2 = cross(n,r1);

% calculate the point coordinates on the circle
xx = r1(1)*cos(phi)+r2(1)*sin(phi);
yy = r1(2)*cos(phi)+r2(2)*sin(phi);
zz = r1(3)*cos(phi)+r2(3)*sin(phi);

x = zeros(2,nPhi);
y = zeros(2,nPhi);
z = zeros(2,nPhi);

x(1,:) = xx + p2(1);
y(1,:) = yy + p2(2);
z(1,:) = zz + p2(3); 
x(2,:) = xx + p1(1);
y(2,:) = yy + p1(2);
z(2,:) = zz + p1(3); 

options = {'FaceColor', clr, 'linestyle', 'none'};
surf(x, y, z, options{:});

end