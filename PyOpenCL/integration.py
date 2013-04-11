import pyopencl as cl
import pyopencl.array as cl_array
import numpy
import numpy.linalg as la
import scipy as sp
from pyopencl.reduction import ReductionKernel
from pyopencl.elementwise import ElementwiseKernel

# Make shure we use the GPU for computations.
platform=cl.get_platforms()
gpu_devices=platform[0].get_devices(device_type=cl.device_type.GPU)
ctx=cl.Context(devices=gpu_devices)
queue = cl.CommandQueue(ctx)

# Integration, one integral, using reduce
length=0.001
vals=cl_array.arange(queue,0,100,length,dtype=numpy.float32)
# Kernel for reduce-code
krnlRed=ReductionKernel(ctx,numpy.float32,neutral="0",
    reduce_expr="a+b",map_expr="get_val(x[i])*%10.3f" % length,
    arguments="__global float *x",
    preamble="""
    float get_val(float x)
    {
        return x*x;
    }
    """)

# Generation of an array where each element is an evaluated integral.
tonum=1000000 # Number of elements.
# Array to send to the GPU.
p_gpu=cl_array.to_device(ctx,queue,sp.linspace(0,tonum,tonum+1).astype(numpy.float32))
res=cl_array.empty_like(p_gpu) # The resultating array

# Elementwise (mapping) kernel.
krnlMap=ElementwiseKernel(ctx,"float *param, float *res", "res[i]=integrate(param[i])",preamble="""
    float integrate(float param)
    {
        float sum=0;
        for (float f=0.0;f<10.0;f+=0.001)
        {
            sum+=(f*f-10*f-param);
        }
        return sum/1000.0;
    }
""")
integrand=krnlRed(vals).get() # Calculate the first integral.
krnlMap(p_gpu,res) # Generate the large array.

# Print results
print "Value of integral: ", integrand
print "x**2 from 0 to 100\n"
print "''Integrand''"
print res.get()