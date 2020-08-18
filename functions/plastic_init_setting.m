function D = plastic_init_setting(mat, ngp, eltype)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Initialize history variables and elastic stiffness matrix
% 
% input:
%       mat:    material definition
%       ngp:    number of integration points
%       eltype: element types,  
%           10 - plane strain element
%           20 - solid element
% 
% SIGMA : Stress for rate-form plasticity
%       : Left Cauchy-Green tensor XB for multiplicative plasticity
% 
% XQ    :   eltype == 10
%               1-4 = Back stress alpha, 5 = Effective plastic strain
%           eltype == 20
%               1-6 = Back stress alpha, 7 = Effective plastic strain
% D     : Elastic stiffness matrix
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

global SIGMA XQ;

E   = mat(2);
nu  = mat(3);
lam = E*nu/(1+nu)/(1-2*nu);
mu  = E/2/(1+nu);

if eltype == 10      % plane element
    
    D = [lam+2*mu    lam        0 ;
         lam         lam+2*mu   0 ; 
         0           0          mu;];
     
    SIGMA = zeros(4,ngp);
    XQ    = zeros(5,ngp);
     
elseif eltype == 20  % solid element

    D = [lam+2*mu lam      lam      0  0  0;
         lam      lam+2*mu lam      0  0  0;
         lam      lam      lam+2*mu 0  0  0;
         0        0        0        mu 0  0;
         0        0        0        0  mu 0;
         0        0        0        0  0  mu];
    
    SIGMA = zeros(6,ngp);
    XQ    = zeros(7,ngp);
end

end