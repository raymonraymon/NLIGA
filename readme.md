# Nonlinear isogeometric analysis - NLIGA
------------------------------------------------------------
This is the README file of the open source MATLAB toolbox: Nonlinear isogeometric analysis (NLIGA). Nonlinear materials including hyperelastic and elastoplastic materials are considered into the code at this stage. Two and three dimensional models are supported. The resulted displacement, stress can be visualized in this toolbox.

------------------------------------------------------------
Please first go to the folder cfiles to run the script compile.m, because some C codes are used in this toolbox to improve the efficiency of the procedure. So a C/C++ compiler is required. 
!!! If you couldn't compile it successfully, don't worry, we can take a detour!

First of all, add the nliga folder to the path. Right-click on the nliga folder and select "Add to Path" -> "Selected Folders and Subfolders". 

To have a quick look of NLIGA, we can directly go to the subfolder 'output' and run the script 'show_examples'. 

The detailed contents of the folders in NLIGA are listed below:
------------------------------------------------------------
1. cfiles
This folder includes the C codes for calculating the two and three dimensional basis functions and derivatives, which will be widely used in other functions. Please run the compile.m file first to generate the mex files which provide an interface between MATLAB and C/C++ codes.

IMPORTANT: If the C codes can't be compiled successfully, please use the functions provided in "nurbs-1.3.12" for calculating derivatives rather than the C codes. The detailed instructions can be found in the function "nurbs_derivatives", which is located in the folder "functions".

------------------------------------------------------------
2. elasticity
A series of benchmark problems in elasticity are considered in this folder. You can find the solving procedures including Poisson problems, plane and solid problems, plate and shell problems (static bending and free vibration). The results obtained in these benchmark problems could be used as references for comparison with results obtained by other methods.

------------------------------------------------------------
3. elastoplasticity
Some elastoplastic examples are presented in this folder. It should be noticed that only small strains in elastoplasticity are supported at this stage. Large rotation and large strains will be supported in the near future. The displacements of control points and the stresses at Gaussian integration points are saved in the corresponding files for visualization. Note that the format of the saved file here is different from that in hyperelasticity. The output file is saved in the folder 'output'.

------------------------------------------------------------
4. /functions
The main functions of NLIGA are given in this folder, including the mesh generation, integration, connectivity, constitutive relations, visualization mesh generation, and plot functions.

------------------------------------------------------------
5. geometries
Several geometries are provided in the folder, where all functions are named with a prefix 'geo' meaning geometrical models. Note that some functions require parameters to control the model. For example, you should provide an initial point and side length to build a square plate.

------------------------------------------------------------
6. hyperelasticity
This folder provides some nonlinear examples in which hyperelastic materials are used. All the scripts are named with a prefix 'ex' meaning example. You can directly run anyone script to find the whole workflow of the NLIGA. The messages including time, time step, iterative step and residual will be shown on the command window, and animation for the results of displacements or stresses will be played in the end. The visualized mesh file '.msh' will be saved in the folder 'output'. If you want to check the animation again, you can run the functions in the folder 'output'.

------------------------------------------------------------
7. multipatch
Considering CAD geometries are usually built with multiple patches, we add some examples to explain how to tackle with multipatch problems. Note that only conforming multi-patch problems are considered and the non-conforming problems will be added in the future by using domain decomposition methods.

------------------------------------------------------------
8. nurbs-1.3.12
This is an excellent opensource toolbox for the creation, manipulation of NURBS. More information can be found on the website, https://octave.sourceforge.io/nurbs/index.html 
The representation of the NURBS and p, k refinement algorithms used in NLIGA are employed from this toolbox.

------------------------------------------------------------
9. output
The mesh files for visualization are saved in this folder. You can review the analytical results by using 'show_examples' to read the .msh files.

------------------------------------------------------------
Features of NLIGA:
1. Nonlinear problems solving by isogeometric analysis; 
2. Hyperelastic and elastoplastic materials are considered;
3. Imposition of the displacement boundary conditons and force boundary conditions;
4. Two and three-dimensional cases are studied;
4. Preprocessing of the geometrical models;
5. Postprocessing of the results includig displacements and stresses;
6. Stress recovery from the Gaussian integration points;
7. Elastic benchmark problems;
8. Multi-patch problems;
9. ...(future update).

------------------------------------------------------------
*** If you have any comments or suggestions, please feel free to contact us! 
- Xiaoxiao Du, 
- Beihang Univeristy,
- duxxhf@gmail.com or 
- duxiaoxiao@buaa.edu.cn. 

*** If you feel it is helpful for your work, please cite our work.




