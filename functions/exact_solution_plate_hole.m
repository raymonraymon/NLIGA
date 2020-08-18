function [exactStress, exactDisplacement] = exact_solution_plate_hole(x, a, E, nu)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Compute exact solution of the infinite plate with a circular hole
%  problem
%  Input:
%    x - coordinates for calculating solutions,[x,y]
%    a - radius of the circular hole
%    E - Youngs modulus
%    nu - Poisson ratio
%  Output:
%    exactStress - calculated exact stress, [S_11, S_22, S_12]
%    exactDisplacement - calculated exact displacement [u_x,u_y]
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

r     = norm(x);
theta = atan(x(2)/x(1));

sigma_xx = 1 - a^2/r^2*( 1.5 *cos(2*theta) + cos(4*theta) ) + 1.5*a^4/r^4 * cos(4*theta);
sigma_yy = -a^2/r^2 *( 0.5* cos(2*theta) - cos(4*theta) ) - 1.5*a^4/r^4 * cos(4*theta);
sigma_xy = -a^2/r^2 *( 0.5* sin(2*theta) + sin(4*theta) ) + 1.5*a^4/r^4 * sin(4*theta);
exactStress = [sigma_xx, sigma_yy, sigma_xy];

mu = E/(2*(1+nu));
kappa = (3-nu)/(1+nu); % for plane stress conditions

u_x = a/(8*mu) * (  r/a * (kappa+1)*cos(theta) + 2*a/r*( (kappa+1)*cos(theta)+cos(3*theta) )- 2*a^3/r^3 *cos(3*theta) );
u_y = a/(8*mu) * (  r/a * (kappa-3)*sin(theta) + 2*a/r*( (1-kappa)*sin(theta)+sin(3*theta) )- 2*a^3/r^3 *sin(3*theta) );
exactDisplacement = [u_x, u_y];

end

