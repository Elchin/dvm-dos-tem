import sys, netCDF4 as nc, numpy as np, glob, os
a,b=sys.argv[1:3]; bad=0; n=0
for st in ['spin000','main']:
  for f in sorted(glob.glob(f'{a}/{st}/output/*.nc')):
    g=f.replace(a,b,1)
    if 'restart' in f or not os.path.exists(g): continue
    A=nc.Dataset(f); B=nc.Dataset(g)
    for v in A.variables:
      x=np.asarray(A[v][:]); y=np.asarray(B[v][:]); n+=1
      if x.shape!=y.shape or not np.array_equal(x,y,equal_nan=True): bad+=1; print('DIFF',st,os.path.basename(f),v)
print('vars compared',n,'differ',bad)
