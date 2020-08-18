function [ pk2, dtan ] = material_hyperelasticity( dim, mat, F )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Calculate constitutive relations  
%  Input:
%    dim - dimension
%    mat - material properties
%    F -  deformation tensor
%  Output:
%    pk2 - Piola-Kirchhoff stress
%    dtan   - tangent stiffness
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

C = F'*F;      % right Cauchy-Green deformation tensor
if dim == 2    % plane element
    C11 = C(1,1); C12 = C(1,2); C21 = C(2,1); C22 = C(2,2);
    C33 = 1 / ( C11*C22 - C12*C21 );   % satisfy det(C) = 1
    I1 = C11+C22+C33;
%     I2 = C11*C22 - C12*C21 + (C11+C22)*C33;
    I1C = [1-C33^2*C22; 1-C33^2*C11; C33^2*C12];
    I2C = [C22+C33-C33^2*(C11+C22)*C22;  C11+C33-C33^2*(C11+C22)*C11; -C12+C33^2*(C11+C22)*C12];

    I1CC = C33^2*[     2*C33*C22^2,  2*C33*C11*C22-1,    -2*C33*C22*C12;
                   2*C33*C11*C22-1,      2*C33*C11^2,    -2*C33*C11*C12;
                    -2*C33*C22*C12,   -2*C33*C11*C12,   2*C33*C12^2+0.5;];

    I2CC = [           2*C33^3*C22^2*(C11+C22)-2*C33^2*C22,   1-2*C33^2*(C11+C22)+2*C33^3*(C11+C22)*C11*C22,            C33^2*C12-2*(C11+C22)*C33^3*C12*C22; ...
             1-2*C33^2*(C11+C22)+2*C33^3*(C11+C22)*C11*C22,             2*C33^3*C11^2*(C11+C22)-2*C33^2*C11,            C33^2*C12-2*(C11+C22)*C33^3*C12*C11; ...               
                       C33^2*C12-2*(C11+C22)*C33^3*C12*C22,             C33^2*C12-2*(C11+C22)*C33^3*C12*C11,  2*(C11+C22)*C33^3*C12^2+C33^2*(C11+C22)/2-0.5;];
    K = mat(2);               
    if mat(1) == 10 && K == 0  % incompressible Neo-Hookean
        A10 = mat(3);
        pk2 = 2*A10*I1C;
        dtan = 4*A10*I1CC;
    elseif mat(1) == 11 && K == 0 % incompressible Mooney-Rivlin
        A10 = mat(3);
        A01 = mat(4);
        pk2 = 2*A10*I1C + 2*A01*I2C;
        dtan = 4*A10*I1CC + 4*A01*I2CC;
    elseif mat(1) == 12 && K == 0 % incompressible Yeoh
        A10 = mat(3); 
        A20 = mat(4);
        A30 = mat(5);
        pk2 = 2*A10*I1C + 4*A20*(I1-3)*I1C + 6*A30*(I1-3)^2*I1C;
        dtan = 4*A10*I1CC + 8*A20*(I1C*I1C') + 8*A20*(I1-3)*I1CC +...
               24*A30*(I1-3)*(I1C*I1C') + 12*A30*(I1-3)^2*I1CC;  
    elseif mat(1) == 13 && K == 0 % incompressible Bidtanerman
        A10 = mat(3); 
        A01 = mat(4);
        A20 = mat(5);
        A30 = mat(6);
        pk2 = 2*A10*I1C + 2*A01*I2C + 4*A20*(I1-3)*I1C + 6*A30*(I1-3)^2*I1C;
        dtan = 4*A10*I1CC + 4*A01*I2CC + 8*A20*(I1C*I1C') + 8*A20*(I1-3)*I1CC +...
            24*A30*(I1-3)*(I1C*I1C') + 12*A30*(I1-3)^2*I1CC; 
    end
elseif dim == 3   % solid element, cosider the nearly-compressible material
    % Right cauchy-green deformation tensor
    C = F'*F;  
    C11 = C(1,1); C22 = C(2,2); C33 = C(3,3); C12 = C(1,2); C23 = C(2,3); C31 = C(3,1);

    % Three invariants of deformation tesnor
    I1 = trace(C);   
    I2 = 0.5*( trace(C)*trace(C)-trace(C*C) );
    I3 = det(C);

    % The third modified invariant
    J3 = sqrt(I3);

    % Voigt representation of first-order derivatives of the three invariants
    I1C = [1;1;1;0;0;0];
    I2C = [C22+C33; C11+C33; C11+C22; -C12; -C23; -C31];
    I3C = [C22*C33-C23*C23; C11*C33-C31*C31; C11*C22-C12*C12; ...
           C23*C31-C12*C33; C12*C31-C11*C23; C12*C23-C22*C31];

    % Matrix representation of second-order derivatives of the three invariants
    I1CC = zeros(6,6);
    I2CC = [0, 1, 1, 0, 0, 0;  1, 0, 1, 0, 0, 0;  1, 1, 0, 0, 0, 0;     ...
            0, 0, 0, -0.5, 0, 0;  0, 0, 0, 0, -0.5, 0;  0, 0, 0, 0, 0, -0.5];
    I3CC = [0, C33, C22, 0, -C23, 0;  C33, 0, C11, 0, 0, -C31;  C22, C11, 0, -C12, 0, 0;  ...
            0, 0, -C12, -C33/2, C31/2, C23/2;  -C23, 0, 0, C31/2, -C11/2, C12/2;  0, -C31, 0, C23/2, C12/2, -C22/2];   

    % First derivatives of the modified invariants   
    J1C = I3^(-1/3)*I1C - 1/3*I1*I3^(-4/3)*I3C;
    J2C = I3^(-2/3)*I2C - 2/3*I2*I3^(-5/3)*I3C;
    J3C = 0.5*I3^(-1/2)*I3C;  

    % Second-order derivatives of the modified invariants
    J1CC = I3^(-1/3)*I1CC + 4/9*I1*I3^(-7/3)*(I3C*I3C') -  1/3*I3^(-4/3)*( I1C*I3C'+I3C*I1C'+I1*I3CC );
    J2CC = I3^(-2/3)*I2CC + 10/9*I2*I3^(-8/3)*(I3C*I3C') - 2/3*I3^(-5/3)*( I2C*I3C'+I3C*I2C'+I2*I3CC );
    J3CC = 0.5*I3^(-0.5)*I3CC - 0.25*I3^(-3/2)*(I3C*I3C');
    K = mat(2);
    if mat(1) == 10  % Neo-Hookean
        % W(J1,J3) = A10*(J1-3)+K/2*(J3-1)^2
        A10 = mat(3);
        pk2 = 2*( A10*J1C + K*(J3-1)*J3C );
        dtan = 4*(A10*J1CC + K*(J3C*J3C') + K*(J3-1)*J3CC);  
    elseif mat(1) == 11 % Mooney-Rivlin
        %  W(J1,J2,J3) = A10*(J1-3)+A01*(J2-3)+K/2*(J3-1)^2
        A10 = mat(3);
        A01 = mat(4);
        pk2 = 2*( A10*J1C + A01*J2C + K*(J3-1)*J3C );
        dtan = 4*(A10*J1CC + A01*J2CC + K*(J3C*J3C') + K*(J3-1)*J3CC);        
    elseif mat(1) == 12 % Yeoh
        % W(J1,J3) = A10*(J1-3)+A20*(J1-3)^2+A30*(J1-3)^3+K/2*(J3-1)^2
        A10 = mat(3);
        A20 = mat(4);
        A30 = mat(5);
        pk2 = 2* ( A10*J1C +2*A20*(J1-3)*J1C+3*A30*(J1-3)^2*J1C+ K*(J3-1)*J3C );
        dtan = 2* ( A10*J1CC +2*A20*J1C*J1C+2*A20*(J1-3)*J1C+6*A30*(J1-3)*J1C*J1C + 3*A30(J1-3)^2*J1CC+ K*(J3C*J3C') + K*(J3-1)*J3CC );
    elseif mat(1) == 13 % Bidtanerman
        % to be updated
    end
end

end

