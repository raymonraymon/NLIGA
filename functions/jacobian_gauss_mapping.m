function jacobian = jacobian_gauss_mapping( elDoma )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Calculate the jacobian of the gauss mapping
%  Input:
%    elDoma - element domain 
%             one-dimensional:   elDoma = [u1, u2]  
%             two-dimensional:   elDoma = [u1, u2, v1, v2] 
%             three-dimensional: elDoma = [u1, u2, v1, v2, w1, w2]  
%  Output:
%    jacobian - jacobian value
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if length(elDoma) == 2        % one dimensional
    jacobian = 0.5 * (elDoma(2)-elDoma(1));
elseif length(elDoma) == 4    % two dimensional
    Jxi = 0.5 * (elDoma(2)-elDoma(1));
    Jeta = 0.5 * (elDoma(4)-elDoma(3));
    jacobian = Jxi * Jeta;
elseif length(elDoma) == 6    % three dimensional
    Jxi = 0.5 * (elDoma(2)-elDoma(1));
    Jeta = 0.5 * (elDoma(4)-elDoma(3));
    Jzeta = 0.5 * (elDoma(6)-elDoma(5));
    jacobian = Jxi * Jeta *Jzeta;
end

end

