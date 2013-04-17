import pyopencl as cl
import pyopencl.array as cl_array
import numpy
import numpy.linalg as la
import scipy as sp
from pyopencl.elementwise import ElementwiseKernel
from pyopencl.reduction import ReductionKernel

size=10

ctx=cl.create_some_context()
queue=cl.CommandQueue(ctx)

a = cl.array.arange(queue, 400, dtype=numpy.float32)
b = cl.array.arange(queue, 400, dtype=numpy.float32)

arr=(numpy.array([i for i in range(size)])).astype(numpy.float32)
gpu_arr=cl_array.to_device(ctx,queue,arr)

krnl = ReductionKernel(ctx, numpy.float32, neutral="0.0f",
        reduce_expr="a+b", map_expr="x[i]*func(y[i])",
        arguments="__global float *x, __global float *y",
        preamble="""
        float func(float x)
            {
                return 1.0f;
            }
        """)

my_dot_prod = krnl(gpu_arr, gpu_arr).get()
print my_dot_prod