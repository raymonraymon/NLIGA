function ga = get_greville_abscissae( geo )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Calculate all the greville abscissae points to related control points
%  Input:
%    geo - geometry structure
%  Output:
%    ga  - coordinates of greville abscissae
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

d = size(geo.number,2);
m = geo.number;
knots = geo.knots;
p = geo.order-1;
if d == 1   % curve
   ga = zeros(m,1);
   for i = 1:m
       ga(i) = sum(knots(i+1:i+p))/p;
   end   
elseif d == 2  % surface
    ga = zeros(m(1)*m(2),2);
    k1 = knots{1,1};
    k2 = knots{1,2};
    count = 1;
    for j = 1:m(2)
        b = sum(k2(j+1:j+p(2)))/p(2);
        for i =1:m(1)
            a = sum(k1(i+1:i+p(1)))/p(1);
            ga(count,:) = [a,b];
            count = count+1;
        end
    end   
elseif d == 3  % volume
    ga = zeros(m(1)*m(2)*m(3),3);
    k1 = knots{1,1};
    k2 = knots{1,2};
    k3 = knots{1,3};
    count = 1;
    for k = 1:m(3)
        c = sum(k3(k+1:k+p(3)))/p(3);
        for j = 1:m(2)
            b = sum(k2(j+1:j+p(2)))/p(2);
            for i = 1:m(1)
                a = sum(k1(i+1:i+p(1)))/p(1);
                ga(count,:) = [a,b,c];
                count = count+1;
            end
        end
    end
    
end

end

