import pyopencl as cl
import pyopencl.array as cl_array
import numpy
import numpy.linalg as la
import scipy as sp
from pyopencl.elementwise import ElementwiseKernel
from pyopencl.reduction import ReductionKernel
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__),'../..'))

from nhqm.calculations import QM as calc

order = 5000
l = 1
j = 1.5
peak_x = 0.17
peak_y = 0.07
k_max = 2.5

contour = calc.gauss_contour([0,3], order)
points, weights = contour
points = points.astype(numpy.float32)
weights = weights.astype(numpy.float32)

ctx=cl.create_some_context()
queue=cl.CommandQueue(ctx)

gpu_points=cl_array.to_device(ctx,queue,points)
gpu_weights=cl_array.to_device(ctx,queue,weights)

krnl = ReductionKernel(ctx, numpy.float32, neutral="0.0f",
        reduce_expr="a+b", map_expr="x[i]*func(y[i])",
        arguments="__global float *x, __global float *y",
        preamble="""
        float func(float x)
            {
                return pow(x,5);
            }
        """)

my_dot_prod = krnl(gpu_weights, gpu_points).get()
print my_dot_prod