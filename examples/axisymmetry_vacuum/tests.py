#!/usr/bin/env python3

from scipy.io import netcdf

f = netcdf.netcdf_file('quasisymmetry_out.axisymmetry_vacuum.nc','r',mmap=False)
R0c = f.variables['R0c'][()]
max_modBinv_sqrt_half_grad_B_colon_grad_B = f.variables['max_modBinv_sqrt_half_grad_B_colon_grad_B'][()]
grad_grad_B_inverse_scale_length = f.variables['grad_grad_B_inverse_scale_length'][()]
f.close()

print('R0c:',R0c)
print('1 / R0c:',1 / R0c)
print('modBinv_sqrt_half_grad_B_colon_grad_B:',max_modBinv_sqrt_half_grad_B_colon_grad_B)
print('grad_grad_B_inverse_scale_length:',grad_grad_B_inverse_scale_length)
assert(abs(1/R0c[0] - max_modBinv_sqrt_half_grad_B_colon_grad_B) < 1.0e-10)
assert(abs(1/R0c[0] - grad_grad_B_inverse_scale_length) < 1.0e-10)
