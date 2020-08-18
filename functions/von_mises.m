function [ vmises ] = von_mises( stress )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Calculate von mises stress from cauchy stress
% stress - Cauchy stress
% vmises - von Mises stress
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%


if size(stress,2) == 3  % two-dimensional
    vmises = sqrt( stress(:,1).^2 - stress(:,1).*stress(:,2) + stress(:,2).^2 + 3*stress(:,3).^2 );
elseif size(stress,2) == 6  % three-dimensional
    vmises = sqrt( 0.5*( (stress(:,1) - stress(:,2)).^2 + (stress(:,2) - stress(:,3)).^2 + (stress(:,3) - stress(:,1)).^2 ) ...
             + 3* ( stress(:,4).^2 + stress(:,5).^2 + stress(:,6).^2 ));
end

end

