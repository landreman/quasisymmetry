#!/usr/bin/env python3

from scipy.io import netcdf

f = netcdf.netcdf_file('quasisymmetry_out.axisymmetry_vacuum.nc','r',mmap=False)
R0c = f.variables['R0c'][()]
grad_grad_B_inverse_scale_length = f.variables['grad_grad_B_inverse_scale_length'][()]
f.close()

print('R0c:',R0c)
print('1 / R0c:',1 / R0c)
print('grad_grad_B_inverse_scale_length:',grad_grad_B_inverse_scale_length)
assert(abs(1/R0c[0] - grad_grad_B_inverse_scale_length) < 1.0e-10)
