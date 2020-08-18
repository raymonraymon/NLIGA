%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  - This is the file to compile C codes into Mex files for MATLAB
%  implementation.
%  - C codes are mainly for calculation of NURBS basis functions and
%  derivatives.
%  - Please make sure you have a C code compiler.
%  - If you try hard but couldn't make it, don't worry, we could take a
%  detour to compute basis functions and their derivatives. Just go to the
%  function 'nurbs_derivatives' in the folder 'functions' to see details.
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

%mex -setup C
mex nrbNurbs1DBasisDerivs.c nrbNurbs.c
mex nrbNurbs1DBasisFunction.c nrbNurbs.c
mex nrbNurbs2DBasisDerivs.c nrbNurbs.c
mex nrbNurbs2DBasisFunction.c nrbNurbs.c
mex nrbNurbs3DBasisDerivs.c nrbNurbs.c
mex nrbNurbs3DBasisFunction.c nrbNurbs.c

