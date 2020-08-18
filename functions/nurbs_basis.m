function R = nurbs_basis( gauPts, geo, mesh )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Calculate nurbs basis dunction
%  Input:
%    gauPts - parameter points, xi or [xi, eta] or [xi, eta, zeta]
%    geo - geometry
%    mesh - mesh structure of the geo
%  Output:
%    R - basis functions
% IMPORTANT: If the C codes can't be compiled successfully, please 
%            follow the instructions below
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

num = length(gauPts);
if num == 1
    R =[];
    % to be added
elseif num == 2
    R =[];
    % to be added
elseif num == 3
    R =[];
    % to be added
end

end

