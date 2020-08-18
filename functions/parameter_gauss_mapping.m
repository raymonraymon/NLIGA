function [ gauPts ] = parameter_gauss_mapping( elDoma, curPts )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Mapping the parameters for current domain to [-1,1]
%  Input:
%    elDoma - element domain 
%             one-dimension:   elDoma = [u1, u2]  
%             two-dimension:   elDoma = [u1, u2, v1, v2] 
%             three-dimension: elDoma = [u1, u2, v1, v2, w1, w2] 
%    curPts -  current coordinates for integration 
%  Output:
%    gauPts - mapped gauss points
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

num = length(curPts);
if isempty(elDoma) || 2*num ~= length(elDoma)
    return;
else
    gauPts = zeros(1,num);
    for i = 1:num
        gauPts(i) = 0.5 * ( ( elDoma(i*2) -elDoma(i*2-1) ) * curPts(i) + ( elDoma(i*2) + elDoma(i*2-1) ) ); 
    end
end

end

