function [ pk2, dtan ] = constitutive_relation( dim, mat, F )
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

if ( mat(1) >= 10 && mat(1) < 20 )            % belongs to hyperelastic materials
    [ pk2, dtan ] = material_hyperelasticity( dim, mat, F );
elseif ( mat(1) >= 20 && mat(1) < 40 )        % belongs to plastic materials
    % to be updated
elseif ( mat(1) >= 40 && mat(1) < 60 )        % belongs to viscous materials
    % to be updated
end

end

