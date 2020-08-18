function [Q, W] = gauss_quadrature(noGpsX, noGpsY, noGpsZ)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Calculate one two and three dimensional gauss integration points and its weights
%  Input:
%    noGpsX - number of gauss points in x direction
%    noGpsY - number of gauss points in y direction
%    noGpsZ - number of gauss points in z direction
%  Output:
%    Q  - coordinates of the integration points
%    w  - weights of the integration points
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

a = -1;    % lower limit of integral interval
b = 1;     % upper limit of integral interval 
if nargin == 1
    dim = 1;
    Q = zeros(noGpsX,1);
    W = zeros(noGpsX,1);
elseif nargin == 2 
    dim = 2;
    Q = zeros(noGpsX*noGpsY,2);
    W = zeros(noGpsX*noGpsY,1);
elseif nargin == 3
    dim = 3;
    Q = zeros(noGpsX*noGpsY*noGpsZ,3);
    W = zeros(noGpsX*noGpsY*noGpsZ,1);
else 
    error('Please input a correct number');
end

if dim == 1    % one dimension
    [Q,W] = gauss_legendre(a,b,noGpsX);
    return;
end

if dim == 2    % two dimension
    [Q1,W1] = gauss_legendre(a,b,noGpsX);
    [Q2,W2] = gauss_legendre(a,b,noGpsY);
    n = 0;
    for j = 1:noGpsY
        for i = 1:noGpsX
            n = n+1;
            Q(n,:) = [Q1(i),Q2(j)];
            W(n) = W1(i)*W2(j);
        end
    end
    return;
end

if dim == 3     % three dimension
    [Q1,W1] = gauss_legendre(a,b,noGpsX);
    [Q2,W2] = gauss_legendre(a,b,noGpsY);
    [Q3,W3] = gauss_legendre(a,b,noGpsZ);
    n = 0;
    for k = 1:noGpsZ
        for j =1:noGpsY
            for i =1:noGpsX
                n = n+1;
                Q(n,:) = [Q1(i),Q2(j),Q3(k)];
                W(n) = W1(i)*W2(j)*W3(k);
            end
        end
    end
    return;
end

end

