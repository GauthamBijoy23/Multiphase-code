IMPORTANT
working on writing cell center value instead of nodes

Can make it faster by storing xai and all variable in geometry only

changed primitive and primitive plus Bc loop now faster

Modified the boundary condtion statement, now there will be no if condtion for BC 

!- if activating YP_P increase allocation of ANU(I, 3) to ANU(I, 4)   

Visocous subroutine in jameson needs to be activaated for the species transport.

Added dellocation of variables 
 
Imp: Re caclulation source qrs is changed fromn up - uq to ug-uq  

KE model is implemented for Jameson 2 phase and 5 phase only

AMUP variation we need for 2 and 5 phase, No need to mix KE with this  


Jamson 18 hllc batten shock tube KE 2 phase GPU 1 

Jameson 16 trying Kepsilon model 2 phase no amup vary

Jameson  started with 5 phase no KE and amup vary, it is not diverging, probably some reason I dont know

Jameson 15 validated for variable amup only 2 phase no KE gpu, not validated is pressure plot is different compared ot old one 
 
Jameson14 it is working with the amup variable, 2phase no KE, not able to validate this case pattern in water phase is diffenrent

Jameson13  Debugging amup 2 phase shock tube check mpi only no gpu and only amup value 
File_name: Jameson_12 Validating Two-phase ke model, amup is constant, only 2 phase actvated, hllcbatten used. Not Working alaone KE with or without AMUP const.  

Made seperate subrouting for initial conditions.

Now implementing KE model.
  
Case9: 5 phase shock tube 0.05% of remainig phase in left side AMUP variable, amup is constant 

Here we are varying amup instead of constant in sourceqrs and  viscougp with ifdefine
akp,amug are now local variable

Cutoff value of the amup is now 0,0001 instead of amug which was very low in some cases

Found one correction in the viscousgp I was using tgf instead of tpf in two places.

One more correction cu and cn loop before the diifu caclulation was wrong corrected it for neq loop

Case : Validating 5 phase activated with amug variable, five phase has cutoff values



3D Multiphase Six-Equation Solver
!================================

This code implements a three-dimensional multiphase six-equation model designed for high-fidelity simulations of two-phase and five-phase flows.
It includes multiple flux computation schemes, GPU acceleration, and MPI-based domain decomposition for large-scale parallel computations.

Flux Schemes
!===========
1 HLLC-ZIEN
2 HLLC-BATTEN
3 HLLC-BATTEN + JAMSON  

These schemes can be selected as per the desired case setup.

Parallelization
!==============

The solver is fully MPI + GPU parallelized.
MPI handles domain-level parallelism, while OpenACC directives enable GPU acceleration for key computational kernels.

How to Run
!=========

1. Set the number of domain decompositions in the file decom.sh.

2. Place the following input files inside the geometry_files directory:
grid.dat → mesh/grid file
bc.in → boundary and execution input file

3. Execute the decomposition script:

./decom.sh

This will decompose the domain into the specified number of parts for parallel execution.

4. Precision and Restart

All CPU-side computations use double precision (16 decimal places).

Updated mpiexcn in gpu. 

Restart files must be written and read in the same precision format to maintain consistency between runs.


OTHER Debugging

5. Boundary Condition Update

The subsonic outflow boundary condition now uses the variable Istart instead of Total for improved accuracy and stability.

Known Issue

A change in pg and pp values has been observed in the first loop of primitive_plus_bc during the third Runge–Kutta (RK3) call when running with GPU acceleration.
The issue is under debugging.

Flux Model Updates

JAMSON scheme has been added for two-phase simulations.

Compile-time flags have been introduced to control two-phase or five-phase compilation using preprocessor definitions (#ifdef FIVE_PHASE).

Changed epsilon value from 10e-8 to 10e-16.


!!!!!Doubt in P+BC we used   do j = 1, neq
    cn(ng,j)=cn(np,j)
   enddo
neq wherever we have used need to recheck 
like in smooth, p, P+BC

!!!!In smooth enerything changed back to 1 to 10 

    


bash run.sh

To complie 
gfortran mesh_to_parts.f90
./a.out

mpmetis mesh_metis.dat 3 -ncommon=3

gfortran mesh_to_parts_check.f90
./a.out

//Metis_check.dat will have infromation of the decomoposed geometry


README.txt
Displaying README.txt.


To compile: 
~~~~~~~~~~~
make 

To clear compilation file:
~~~~~~~~~~~~~~~~~~~~~~~~~
make clean


About:
~~~~~

Multiphase 3D code with ACC parallel

Flux Schmes:
1 HLLC ZIEN
2 HLLC BATTEN

Will write multiple output and restart files for "niter" times


Serial
~~~~~~

In make file remove -acc flag 


Gprof
~~~~~
Compile code with gprof add -pg during the compilation
run the executable and 
gprof a.out > outpt.txt
 

Profiling eith NSIGHT
~~~~~~~~~~~~~~~~~~~~~

nsys profile ./a.out
open .rep file in nsight

Update in the code
~~~~~~~~~~~~~~~~~~
1. Uncommented some parts in hllc batten need to repeat in hllc zien also. Check it 
    nt = 0
    if (neigh > neles) nt = ntype(neigh - neles)
    
2. Changed sequence of viscousgp call, it should be called just after hllcbatten       
