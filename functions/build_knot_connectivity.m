function [ elNodeCnt, elDoma ] = build_knot_connectivity( knot )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Build node element contivity from one dimensional knot vector
%  Input:
%    knots - knot vector
%  Output:
%    elNodeCnt - node connectivity
%    elDoma - element parameter domain
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

uniqueKnot = unique(knot);      % unique knots
nelem = length(uniqueKnot)-1;    % number of elements
p = length(find(knot == knot(1)))-1;  % degree
elNodeCnt = zeros(nelem, p+1);
elDoma = zeros(nelem, 2);
temp = 1;
for i = 1:nelem
    elDoma(i,1) = uniqueKnot(i);
    elDoma(i,2) = uniqueKnot(i+1);
    if i == 1  
        elNodeCnt(i,:) = 1:p+1;
    else
        repeatKnotNum = length(find(knot == uniqueKnot(i)));
        temp = temp + repeatKnotNum;
        elNodeCnt(i,:) = temp:temp+p;
    end
end

end

