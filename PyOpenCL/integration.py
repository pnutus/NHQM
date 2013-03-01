import pyopencl as cl
import pyopencl.array as cl_array
import numpy
import numpy.linalg as la
import scipy as sp
from pyopencl.reduction import ReductionKernel
from pyopencl.elementwise import ElementwiseKernel

# Must-have-variables for OpenCL
ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)

length=0.001
vals=cl_array.arange(queue,0,100,length,dtype=numpy.float32)

krnlRed=ReductionKernel(ctx,numpy.float32,neutral="0",
    reduce_expr="a+b",map_expr="get_val(x[i])*%10.3f" % length,
    arguments="__global float *x",
    preamble="""
    float get_val(float x)
    {
        return x*x;
    }
    """)

tonum=1000000
p_gpu=cl_array.to_device(ctx,queue,sp.linspace(0,tonum,tonum+1).astype(numpy.float32))
res=cl_array.empty_like(p_gpu)
    
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
integrand=krnlRed(vals).get()
krnlMap(p_gpu,res)

print "Value of integral: ", integrand
print "x**2 from 0 to 100\n"
print "''Integrand''"
print res.get()