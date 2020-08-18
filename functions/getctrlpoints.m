function points = getctrlpoints(nurbs)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% get the coordinates of control points of nurbs curve, surface or volume
% 
% acceptable forms:
%       points = getctrlpoints(nurbs)
% 
% inputs:
%       nurbs: nurbs curve, surface or volume
% 
% outputs:
%       points: coordinates of control points, 
%       for curve,
%           points(3, nurbs.number) = zeros(3,nurbs.number(1),nurbs.number(2),nurbs.number(3));
%           where points(i, :) = [xi;yi;zi]
%       for surface,
%           points(3, nurbs.number(1),nurbs.number(2)) = zeros(3,nurbs.number(1),nurbs.number(2),nurbs.number(3));
%           where points(i, :, :) = [xi;yi;zi]
%       for volume,
%           points(3,nurbs.number(1),nurbs.number(2),nurbs.number(3)) = zeros(3,nurbs.number(1),nurbs.number(2),nurbs.number(3));
%           where points(i, :, :, :) = [xi;yi;zi]
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
tic;

if (iscell (nurbs.knots))                                       % nurbs volume or surface
    if (size (nurbs.knots,2) == 2)                              % nurbs surface
        
        points = zeros(3,nurbs.number(1),nurbs.number(2));
        c = 1;
        for j = 1:nurbs.number(2)
            for i = 1:nurbs.number(1)
                points(:,i,j) = nurbs.coefs(1:3,i,j)/nurbs.coefs(4,i,j);
                c = c+1;
            end
        end
        
    elseif (size (nurbs.knots,2) == 3)                          % nurbs volume
        
        points = zeros(3,nurbs.number(1),nurbs.number(2),nurbs.number(3));
        c = 1;
        for k = 1:nurbs.number(3)
            for j = 1:nurbs.number(2)
                for i = 1:nurbs.number(1)
                    points(:,i,j,k) = nurbs.coefs(1:3,i,j,k)/nurbs.coefs(4,i,j,k);
                    c = c+1;
                end
            end
        end  
        
    end
else                                                            % nurbs curve
    
    points = zeros(3,nurbs.number);
    for i = 1:nurbs.number
        points(:,i) = nurbs.coefs(1:3,i)/nurbs.coefs(4,i);
    end
    
end

% disp(['get control points time ', num2str(toc)]);

end