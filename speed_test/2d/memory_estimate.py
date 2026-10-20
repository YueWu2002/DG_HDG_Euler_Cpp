#!/usr/bin/env python3
"""Major numerical pools; RSS also includes allocator/library bookkeeping."""
import argparse
import json
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--cells',type=int,default=256);p.add_argument('--k',type=int,default=8)
a=p.parse_args();n=a.k+1;C=a.cells*a.cells;m=4*n*n;f=4*n;t=4*f
GiB=1024**3
pools=dict(dg_dense_lu=C*m*m*8,
    hdg_recovery_X=C*m*t*8,hdg_local_schur=C*t*t*8,hdg_recovery_y=C*m*8,
    hdg_face_jacobi_lu=2*C*f*f*8,dg_gmres_basis=C*m*21*8,hdg_gmres_basis=2*C*f*21*8,
    volume_coefficient_cache=C*n*n*3*8,face_coefficient_cache=2*C*n*12*8)
print(json.dumps({name:round(value/GiB,4) for name,value in pools.items()},indent=2))
