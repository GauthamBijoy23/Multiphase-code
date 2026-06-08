#to change no. of partitions
# 2 changes here  (mpmetis)
# 2 changes in mesh_to_parts_check ( metis_npart and epart)
# 2 changes in preproc file (npartition and reading metis file)
#--------------------------------------------
#!/bin/bash

NPARTS=1   # set number of partitions only once

echo "Number of partitions = $NPARTS"

# Compile and run mesh preprocessing
gfortran metis_files/mesh_to_parts.f90 -o mesh_to_parts
./mesh_to_parts

rm -f metis_files/mesh_metis.dat.epart.* metis_files/mesh_metis.dat.npart.*

if [ "$NPARTS" -eq 1 ]; then
    echo "Single partition detected: skipping METIS"

    # Ensure output directories exist
    mkdir -p geometry_files

    # Copy geometry and BC files
    cp geometry_files/bc.in   geometry_files/bc_proc0000.in
    cp geometry_files/grid.dat geometry_files/grid_proc0000.dat

    # Create proc_bc0000.in with required content
    echo "0    0" > geometry_files/proc_bc0000.in

else
    echo "Running METIS partitioning"
    mpmetis metis_files/mesh_metis.dat $NPARTS -ncommon=3

    gfortran metis_files/mesh_to_parts_check.f90 -o mesh_to_parts_check
    ./mesh_to_parts_check

    rm -f geometry_files/proc_bc* geometry_files/grid_proc* geometry_files/bc_proc*

    gfortran preProc_file/f23.f90 -fbacktrace -o preProc
    ./preProc
fi

# Cleanup
rm -f mesh_to_parts mesh_to_parts_check preProc metis_var.mod


#g++ -std=c++17 -o find_cell preProc_file/Probe_location_track.C




