function [stress, displacement]= exact_stress_curved_beam( x, a, b, nu, E )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Compute exact stress for the curved beam example given in 'ex_elast_curved_beam'
%  Input:
%    x - coordinates in form of [x,y]
%    a - inner radius
%    b - outer radius
%    nu - Poisson ratio
%    E - Young's modulus
%  Output:
%    stress -  exact stress
%    displacement - exact displacement
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

N = a^2-b^2 + (a^2+b^2)*log(b)/a;
r = sqrt(x(1)^2+x(2)^2);
theta =atan(x(2)/x(1));
pstress = zeros(3,1);
pstress(1) = 1/N * ( r + a^2*b^2/r^3 - (a^2+b^2)/r )*sin(theta);  % exact polar stress
pstress(2) = 1/N * ( 3*r - a^2*b^2/r^3 - (a^2+b^2)/r )*sin(theta);
pstress(3) = -1/N * ( r + a^2*b^2/r^3 - (a^2+b^2)/r )*cos(theta);

% construct transformation matrix for converting stress in polar
% coordinates to normal coordinates
T = [cos(theta)^2, sin(theta)^2, -sin(2*theta);
     sin(theta)^2, cos(theta)^2,  sin(2*theta);
     sin(2*theta)/2, -sin(2*theta)/2, cos(2*theta);];
stress =  T*pstress;
 
K = 1/(N*E) * ( (1-3*nu)*a^2/2 - b^2*(1+nu)/2 - (a^2+b^2)*(1-nu)*log(a) );

ur = 1/(N*E) * ( ((1-3*nu)*r^2/2 - a^2*b^2*(1+nu)/(2*r^2) - (a^2+b^2)*(1-nu)*log(r))*sin(theta)...
    + (a^2+b^2)*(2*theta-pi)*cos(theta) ) - K*sin(theta);
utheta = -1/(N*E) * ( ((5+nu)*r^2/2 - a^2*b^2*(1+nu)/(2*r^2) + (a^2+b^2)*((1-nu)*log(r) ...
    + (1+nu)))*cos(theta) + (a^2+b^2)*(2*theta-pi)*sin(theta) ) - K*cos(theta);

ux = ur*cos(theta) - utheta*sin(theta);
uy = ur*sin(theta) + utheta*cos(theta);
displacement = [ux; uy];

end

