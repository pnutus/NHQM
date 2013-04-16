import pyopencl as cl
import pyopencl.array as cl_array
import numpy
import numpy.linalg as la
import scipy as sp
from pyopencl.elementwise import ElementwiseKernel


size=10

ctx=cl.create_some_context()
queue=cl.CommandQueue(ctx)

arr=(numpy.array([i for i in range(size)])).astype(numpy.float32)
gpu_arr=cl_array.to_device(ctx,queue,arr)
gpu_res=cl_array.empty_like(gpu_arr)

kernel=ElementwiseKernel(ctx,"float *num, float *array, float *res", "res[i]=func(num[i], array[1])",\
                            preamble="""
                            float func(float a, float b)
                            {
                                return a+*(&b+(int)a);
                            }
                            """)
kernel(gpu_arr,gpu_arr,gpu_res)
print gpu_res.get()


# Does not work currently...
# It should return 0,2,4,6, etc, it doesn't