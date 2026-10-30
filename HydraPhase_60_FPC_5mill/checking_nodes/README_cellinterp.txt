FILES NEEDED: 

1. gridfile = 'grid_proc0000.dat'
2. prefile  = 'qcell_preinterp_0000.dat'
3. mapfile  = 'donorinfo_cell.dat'
4. postfile = 'qcell_postinterp_0000.dat' 

CHECK FLOW : 

1. checkcelldon.f90 - obj: checkcelldon - checks if receptor inside donor boundaries - lambda(1to4) = dist from barycenter
                                        - Outputs if any rec not inside donor and max lambda
                                        - checks actual interp value diff (should be 0, bcs straight copy from donr. to recp.
