function V = voigt( M )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Transform second order matrix form to voigt notation
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if size(M,1) == 2
    V = [M(1,1) M(2,2) M(1,2)]';
elseif size(M,1) == 3
    V = [M(1,1) M(2,2) M(3,3) M(1,2) M(2,3) M(1,3)]';
end

end

