Before using the M files, you should first compile some C codes implemented as MEX files.
You have two ways to compile:

	1. Open Matlab and in the command line, type the following commands:
	>>compile
	2. Open Matlab and enter the file named "compile.m", click the run button

Actually compile.m contains commands to compile MEX files.
Our MEX files were compiled on win7 with Microsoft Visual C++ 2010 (C) and Ubuntu 18.04(Linux). We have not tested with other compilers.

These files can be found in folder "cfiles".
1. nrbNurbs.h and nrbNurbs.c
	Declare some functions can be called by other C files and
	they do not generate corresponding MEX files after compiling
2. nrbNurbs1DBasisDerivs.c
	Return the 1D NURBS basis functions and first derivatives to matlab
3. nrbNurbs1DBasisFunction.c
	Return the 1D NURBS basis functions to matlab
4. nrbNurbs2DBasisDerivs.c
	Return the 2D NURBS basis functions and first derivatives to matlab
5. nrbNurbs2DBasisFunction.c
	Return the 2D NURBS basis functions to matlab
6. nrbNurbs3DBasisDerivs.c
	Return the 3D NURBS basis functions and first derivatives to matlab
7. nrbNurbs3DBasisFunction.c
	Return the 3D NURBS basis functions to matlab
