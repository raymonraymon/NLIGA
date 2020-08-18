function Q = local2global_voigt( q )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  convert the local-to-global matrix in tensor-form to voigt form
%  Input:
%    q - local-to-gobal transform matrix in tensor form, 3x3 matrix
%  Output:
%    Q  - local-to-gobal transform matrix in voigt form, 6x6 matrix
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% six component in the first row
Q(1,1) = q(1,1) * q(1,1); Q(1,2) = q(1,2) * q(1,2); Q(1,3) = q(1,3) * q(1,3);
Q(1,4) = 2 * q(1,1) * q(1,2); Q(1,5) = 2 * q(1,2) * q(1,3); Q(1,6) = 2 * q(1,1) * q(1,3);

% six component in the second row
Q(2,1) = q(2,1) * q(2,1); Q(2,2) = q(2,2) * q(2,2); Q(2,3) = q(2,3) * q(2,3);
Q(2,4) = 2 * q(2,1) * q(2,2); Q(2,5) = 2 * q(2,2) * q(2,3); Q(2,6) = 2 * q(2,1) * q(2,3);

% six component in the third row
Q(3,1) = q(3,1) * q(3,1); Q(3,2) = q(3,2) * q(3,2); Q(3,3) = q(3,3) * q(3,3);
Q(3,4) = 2 * q(3,1) * q(3,2); Q(3,5) = 2 * q(3,2) * q(3,3); Q(3,6) = 2 * q(3,1) * q(3,3);

% six component in the fourth row
Q(4,1) = q(1,1) * q(2,1); Q(4,2) = q(1,2) * q(2,2); Q(4,3) = q(1,3) * q(2,3);
Q(4,4) = q(1,2) * q(2,1) + q(1,1) * q(2,2); Q(4,5) = q(1,3) * q(2,2) + q(1,2) * q(2,3); 
Q(4,6) = q(1,3) * q(2,1) + q(1,1) * q(2,3);

% six component in the fifth row
Q(5,1) = q(2,1) * q(3,1); Q(5,2) = q(2,2) * q(3,2); Q(5,3) = q(2,3) * q(3,3);
Q(5,4) = q(2,2) * q(3,1) + q(2,1) * q(3,2); Q(5,5) = q(2,3) * q(3,2) + q(2,2) * q(3,3); 
Q(5,6) = q(2,3) * q(3,1) + q(2,1) * q(3,3);

% six component in the sixth row
Q(6,1) = q(1,1) * q(3,1); Q(6,2) = q(1,2) * q(3,2); Q(6,3) = q(1,3) * q(3,3);
Q(6,4) = q(1,2) * q(3,1) + q(1,1) * q(3,2); Q(6,5) = q(1,3) * q(3,2) + q(1,2) * q(3,3); 
Q(6,6) = q(1,3) * q(3,1) + q(1,1) * q(3,3);



% % six component in the first row
% Q(1,1) = q(1,1) * q(1,1); Q(1,2) = q(1,2) * q(1,2); Q(1,3) = q(1,3) * q(1,3);
% Q(1,4) = q(1,1) * q(1,2); Q(1,5) = q(1,2) * q(1,3); Q(1,6) = q(1,1) * q(1,3);
% 
% % six component in the second row
% Q(2,1) = q(2,1) * q(2,1); Q(2,2) = q(2,2) * q(2,2); Q(2,3) = q(2,3) * q(2,3);
% Q(2,4) = q(2,1) * q(2,2); Q(2,5) = q(2,2) * q(2,3); Q(2,6) = q(2,1) * q(2,3);
% 
% % six component in the third row
% Q(3,1) = q(3,1) * q(3,1); Q(3,2) = q(3,2) * q(3,2); Q(3,3) = q(3,3) * q(3,3);
% Q(3,4) = q(3,1) * q(3,2); Q(3,5) = q(3,2) * q(3,3); Q(3,6) = q(3,1) * q(3,3);
% 
% % six component in the fourth row
% Q(4,1) = 2 * q(1,1) * q(2,1); Q(4,2) = 2 * q(1,2) * q(2,2); Q(4,3) = 2 * q(1,3) * q(2,3);
% Q(4,4) = q(1,2) * q(2,1) + q(1,1) * q(2,2); Q(4,5) = q(1,3) * q(2,2) + q(1,2) * q(2,3); 
% Q(4,6) = q(1,3) * q(2,1) + q(1,1) * q(2,3);
% 
% % six component in the fifth row
% Q(5,1) = 2 * q(2,1) * q(3,1); Q(5,2) = 2 * q(2,2) * q(3,2); Q(5,3) = 2 * q(2,3) * q(3,3);
% Q(5,4) = q(2,2) * q(3,1) + q(2,1) * q(3,2); Q(5,5) = q(2,3) * q(3,2) + q(2,2) * q(3,3); 
% Q(5,6) = q(2,3) * q(3,1) + q(2,1) * q(3,3);
% 
% % six component in the sixth row
% Q(6,1) = 2 * q(1,1) * q(3,1); Q(6,2) = 2 * q(1,2) * q(3,2); Q(6,3) = 2 * q(1,3) * q(3,3);
% Q(6,4) = q(1,2) * q(3,1) + q(1,1) * q(3,2); Q(6,5) = q(1,3) * q(3,2) + q(1,2) * q(3,3); 
% Q(6,6) = q(1,3) * q(3,1) + q(1,1) * q(3,3);

end

