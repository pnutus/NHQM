import pyopencl as cl
import pyopencl.array as cl_array
import numpy
import numpy.linalg as la
import scipy as sp
from pyopencl.elementwise import ElementwiseKernel
import time

# Make shure we use the GPU for computations.
platform=cl.get_platforms()
gpu_devices=platform[0].get_devices(device_type=cl.device_type.GPU)
ctx=cl.Context(devices=gpu_devices) # Context
queue = cl.CommandQueue(ctx) # Queue

# Timing, to see how fast the code is now.
t1=time.time()

# Create all the necessary arrays.
side_length=1000
x_vector=numpy.array([i%side_length for i in range(side_length**2)])
y_vector=numpy.array([(i-i%side_length)/side_length for i in range(side_length**2)])
gpu_matrix_x=cl_array.to_device(ctx,queue,(x_vector).astype(numpy.float32))
gpu_matrix_y=cl_array.to_device(ctx,queue,(y_vector).astype(numpy.float32))
gpu_matrix_res=cl_array.empty_like(gpu_matrix_x)

# Kernel to generate an identity matrix.
krnl_identity_matrix=ElementwiseKernel(ctx,"float *x, float *y, float *res",
    "res[i]=get_element(x[i],y[i])",preamble="""
    float get_element(float x, float y)
    {
        if ((x-y)<0)
        {
            x*=-1.0;
            y*=-1.0;
        }
        return (x-y)<0.0001;
    }
    """)
 
# Kernel to generate a matrix through integrals depending on row and column. 
# The integrand is given in the C-function f().
krnl_gaussian_matrix=ElementwiseKernel(ctx,"float *x, float *y, float start, float end, float *res",
    "res[i]=get_element(x[i],y[i],start,end)",preamble="""
    float f(float r,float x,float y)
    {
        return r;
    }
    float get_element(float x, float y, float start, float end)
    {
        float w[10]={0.2955242247147529,0.2955242247147529,0.2692667193099963,0.2692667193099963,
            0.2190863625159820,	0.2190863625159820,	0.1494513491505806,	0.1494513491505806,
            0.0666713443086881,	0.0666713443086881};
        float xx[10]={-0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472,
            -0.6794095682990244,0.6794095682990244,-0.8650633666889845,0.8650633666889845,
            -0.9739065285171717,0.9739065285171717};
        float sum=0.0;
        for (int i=0;i<10;i++)
        {
            sum += w[i]*f(((end-start)*(xx[i]+1.0)/2.0 + start),x,y);
        }
        return (end-start)*sum/2.0;
    }
    """)
    
# Chose which kernel to execute.
#krnl_identity_matrix(gpu_matrix_x,gpu_matrix_y,gpu_matrix_res)
krnl_gaussian_matrix(gpu_matrix_x,gpu_matrix_y,0,10,gpu_matrix_res)
print gpu_matrix_res.get()

# Timing, again.
t2=time.time()
print "Time: ", (t2-t1)