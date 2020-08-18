function rotor = geo_micro_motor_rotor_8patch( r1, r2, r3, alpha2, alpha3 )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Using eight patches to build one-quarter of a micro-motor rotor
% Input:
%   r1, r2, r3 - radius 
%   alpha2,alpha3 - angles 
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 5    % default parameters
    r1 = 8*10^(-6);
    r2 = 25*10^-6;
    r3 = 50*10^-6;
    rad = pi/180;
    alpha1 = 10*rad;
    alpha2 = 20*rad;
    alpha3 = 25*rad;
else
    rad = pi/180;
    alpha1 = alpha2/2*rad;
    alpha2 = alpha2*rad;
    alpha3 = alpha3*rad;
end

w1 = cos(alpha1/2);
w2 = cos(alpha2/2);
w3 = cos(alpha3/2);
% Build patch 1
coefs1 = zeros(4,3,2);
coefs1(:,:,1) = [0, r1,0,1; r1*tan(alpha1/2)*w1, r1*w1, 0, w1; r1*sin(alpha1), r1*cos(alpha1),0,1]';
coefs1(:,:,2) = [0, r2,0,1; r2*tan(alpha1/2)*w1, r2*w1, 0, w1; r2*sin(alpha1), r2*cos(alpha1),0,1]';
knots1{1} = [0 0 0 1 1 1];
knots1{2} = [0 0 1 1];
rotor{1,1} = nrbmak(coefs1, knots1);

% Build patch 2
coefs2 = zeros(4,3,2);
coefs2(:,:,1) = [ r1*sin(alpha1), r1*cos(alpha1),0,1; 
                  r1/cos(alpha3/2)*sin(alpha1+alpha3/2)*w3, r1/cos(alpha3/2)*cos(alpha1+alpha3/2)*w3,0,w3; 
                  r1*sin(alpha1+alpha3), r1*cos(alpha1+alpha3),0,1]';
coefs2(:,:,2) = [ r2*sin(alpha1), r2*cos(alpha1),0,1; 
                  r2/cos(alpha3/2)*sin(alpha1+alpha3/2)*w3, r2/cos(alpha3/2)*cos(alpha1+alpha3/2)*w3,0,w3; 
                  r2*sin(alpha1+alpha3), r2*cos(alpha1+alpha3),0,1]';
knots2{1} = [0 0 0 1 1 1];
knots2{2} = [0 0 1 1];
rotor{1,2} = nrbmak(coefs2, knots2);

% Build patch 3
coefs3 = zeros(4,3,2);
coefs3(:,:,1) = [ r1*sin(alpha1+alpha3), r1*cos(alpha1+alpha3),0,1; 
                  r1/cos(alpha2/2)*sin(alpha1+alpha3+alpha2/2)*w2, r1/cos(alpha2/2)*cos(alpha1+alpha3+alpha2/2)*w2,0,w2; 
                  r1*cos(alpha1+alpha3), r1*sin(alpha1+alpha3),0,1; ]';
coefs3(:,:,2) = [ r2*sin(alpha1+alpha3), r2*cos(alpha1+alpha3),0,1; 
                  r2/cos(alpha2/2)*sin(alpha1+alpha3+alpha2/2)*w2, r2/cos(alpha2/2)*cos(alpha1+alpha3+alpha2/2)*w2,0,w2; 
                  r2*cos(alpha1+alpha3), r2*sin(alpha1+alpha3),0,1; ]';
knots3{1} = [0 0 0 1 1 1];
knots3{2} = [0 0 1 1];
rotor{1,3} = nrbmak(coefs3, knots3);

% Build patch 4
coefs4 = zeros(4,3,2);
coefs4(:,:,1) = [ r1*cos(alpha1+alpha3), r1*sin(alpha1+alpha3),0,1; 
                  r1/cos(alpha3/2)*cos(alpha1+alpha3/2)*w3, r1/cos(alpha3/2)*sin(alpha1+alpha3/2)*w3,0,w3; 
                  r1*cos(alpha1), r1*sin(alpha1),0,1]';
coefs4(:,:,2) = [ r2*cos(alpha1+alpha3), r2*sin(alpha1+alpha3),0,1; 
                  r2/cos(alpha3/2)*cos(alpha1+alpha3/2)*w3, r2/cos(alpha3/2)*sin(alpha1+alpha3/2)*w3,0,w3; 
                  r2*cos(alpha1), r2*sin(alpha1),0,1]';
knots4{1} = [0 0 0 1 1 1];
knots4{2} = [0 0 1 1];
rotor{1,4} = nrbmak(coefs4, knots4);

% Build patch 5
coefs5 = zeros(4,3,2);
coefs5(:,:,1) = [ r1*cos(alpha1), r1*sin(alpha1),0,1; 
                  r1*w1, r1*tan(alpha1/2)*w1,0,w1; 
                  r1, 0, 0, 1]';
coefs5(:,:,2) = [ r2*cos(alpha1), r2*sin(alpha1),0,1; 
                  r2*w1, r2*tan(alpha1/2)*w1,0,w1; 
                  r2, 0, 0, 1]';
knots5{1} = [0 0 0 1 1 1];
knots5{2} = [0 0 1 1];
rotor{1,5} = nrbmak(coefs5, knots5);

% Build patch 6
coefs6 = zeros(4,3,2);
coefs6(:,:,1) = [0, r2,0,1; r2*tan(alpha1/2)*w1, r2*w1, 0, w1; r2*sin(alpha1), r2*cos(alpha1),0,1]';
coefs6(:,:,2) = [0, r3,0,1; r3*tan(alpha1/2)*w1, r3*w1, 0, w1; r3*sin(alpha1), r3*cos(alpha1),0,1]';
knots6{1} = [0 0 0 1 1 1];
knots6{2} = [0 0 1 1];
rotor{1,6} = nrbmak(coefs6, knots6);

% Build patch 7
coefs7 = zeros(4,3,2);
coefs7(:,:,1) = [ r2*sin(alpha1+alpha3), r2*cos(alpha1+alpha3),0,1; 
                  r2/cos(alpha2/2)*sin(alpha1+alpha3+alpha2/2)*w2, r2/cos(alpha2/2)*cos(alpha1+alpha3+alpha2/2)*w2,0,w2; 
                  r2*cos(alpha1+alpha3), r2*sin(alpha1+alpha3),0,1; ]';
coefs7(:,:,2) = [ r3*sin(alpha1+alpha3), r3*cos(alpha1+alpha3),0,1; 
                  r3/cos(alpha2/2)*sin(alpha1+alpha3+alpha2/2)*w2, r3/cos(alpha2/2)*cos(alpha1+alpha3+alpha2/2)*w2,0,w2; 
                  r3*cos(alpha1+alpha3), r3*sin(alpha1+alpha3),0,1; ]';
knots7{1} = [0 0 0 1 1 1];
knots7{2} = [0 0 1 1];
rotor{1,7} = nrbmak(coefs7, knots7);

% Build patch 8
coefs8 = zeros(4,3,2);
coefs8(:,:,1) = [ r2*cos(alpha1), r2*sin(alpha1),0,1; 
                  r2*w1, r2*tan(alpha1/2)*w1,0,w1; 
                  r2, 0, 0, 1]';
coefs8(:,:,2) = [ r3*cos(alpha1), r3*sin(alpha1),0,1; 
                  r3*w1, r3*tan(alpha1/2)*w1,0,w1; 
                  r3, 0, 0, 1]';
knots8{1} = [0 0 0 1 1 1];
knots8{2} = [0 0 1 1];
rotor{1,8} = nrbmak(coefs8, knots8);

% degeree elevate  of patch 1 - 8
for i = 1:8
    rotor{1,i} = nrbdegelev(rotor{1,i},[0,1]);
end

% insert knots
RefinementX = 3;    % the number of knots inseted in u direction 
RefinementY = 3*RefinementX;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
rotor{1,1} = nrbkntins(rotor{1,1}, {iuknots ivknots});
rotor{1,5} = nrbkntins(rotor{1,5}, {iuknots ivknots});
rotor{1,6} = nrbkntins(rotor{1,6}, {iuknots ivknots});
rotor{1,8} = nrbkntins(rotor{1,8}, {iuknots ivknots});

RefinementX2 = 2*RefinementX;    % the number of knots inseted in u direction 
RefinementY2 = RefinementY;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX2+1):1/(RefinementX2+1):RefinementX2/(RefinementX2+1);
ivknots = 1/(RefinementY2+1):1/(RefinementY2+1):RefinementY2/(RefinementY2+1);
rotor{1,2} = nrbkntins(rotor{1,2}, {iuknots ivknots});
rotor{1,3} = nrbkntins(rotor{1,3}, {iuknots ivknots});
rotor{1,4} = nrbkntins(rotor{1,4}, {iuknots ivknots});
rotor{1,7} = nrbkntins(rotor{1,7}, {iuknots ivknots});

figure
for i = 1:length(rotor)
    hold on;
    plot_nurbs(rotor{1,i}, 0,1);
end
axis([-inf,inf,-inf,inf,-inf,inf]);
view(2);

end

