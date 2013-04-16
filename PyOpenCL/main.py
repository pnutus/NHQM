import pyopencl as cl
import pyopencl.array as cl_array
import numpy
import numpy.linalg as la
import scipy as sp
from pyopencl.elementwise import ElementwiseKernel
from gen_matrix import GenMatrix

size=1000

gm=GenMatrix()
gm.load_potential("woods-saxon.cl")
gm.set_method("mom_space_real.cl")
gm.allocate_space(size,numpy.float32)
gm.combine_kernel("1.0")
gm.execute_kernel(0,6.0)
print gm.get_results().reshape((size,size))


# [eigs,_]=la.eig(gm.get_results().reshape((size,size)))
# indexes = eigs.argsort()
# eigs = sp.real_if_close(eigs[indexes])
# print eigs[0]