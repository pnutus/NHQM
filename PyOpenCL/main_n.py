import pyopencl as cl
import pyopencl.array as cl_array
import numpy
import numpy.linalg as la
import scipy as sp
from pyopencl.elementwise import ElementwiseKernel
from gen_matrix_n import GenMatrix
import sys,os.path
sys.path.append(os.path.join(os.path.dirname(__file__),'..'))
from nhqm.calculations import QM as calc

size=1000

gm=GenMatrix()
gm.load_potential("woods-saxon_complex.cl")
gm.set_method("mom_space_complex.cl")
gm.allocate_space(numpy.complex64,calc.triangle_contour(0.17,0.07,2.5,20))
gm.combine_kernel("1.0")
gm.execute_kernel(0,6.0)
# print gm.get_results().reshape((size,size))
print gm.get_results()


# [eigs,_]=la.eig(gm.get_results().reshape((size,size)))
# indexes = eigs.argsort()
# eigs = sp.real_if_close(eigs[indexes])
# print eigs[0]