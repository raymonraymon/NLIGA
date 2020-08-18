function ccy = pk2cauchy( pk2, f )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Transform second Piola-Kirchhoff stress to cauchy stress
%  formula: ccy = fsf'/j
%  Input:
%    pk2 - second Piola-Kirchhoff stress
%    f -  deformation tensor
%  Output:
%    ccy - Cauchy stress
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if length(pk2) == 6       % three dimensional
    pk2 = [pk2(1), pk2(4), pk2(6);
           pk2(4), pk2(2), pk2(5);
           pk2(6), pk2(5), pk2(3)];
elseif length(pk2) == 3   % two dimensional
    pk2 = [pk2(1), pk2(3);
           pk2(3), pk2(2)];
end

% cauchy = f*pk2*f'/det(f);
ccy = f*pk2*f';

ccy = voigt(ccy);

end
