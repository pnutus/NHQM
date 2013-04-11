###################
# Comments
###################

import pyopencl as cl
import pyopencl.array as cl_array
import numpy
import numpy.linalg as la
import scipy as sp
from pyopencl.elementwise import ElementwiseKernel

class GenMatrix:
    # Set up OpenCL to run at the GPU.
    def __init__(self):
        platform=cl.get_platforms()
        gpu_devices=platform[0].get_devices(device_type=cl.device_type.GPU)
        self.ctx=cl.Context(devices=gpu_devices)
        self.queue = cl.CommandQueue(self.ctx)
    # Load the desired potential, specified in filename.
    def load_potential(self, filename="simple.cl"):
        self.potential="".join(open(filename,'r').readlines())
    # Choose which calculation method to use.
    def set_method(self,method="mom_space_real"):
        self.method="".join(open(method+".cl",'r').readlines())
    def allocate_space(self,size):
        self.size=size
        host_matrix=(numpy.array([i for i in range(size**2)])).astype(numpy.int32)
        self.gpu_matrix=cl_array.to_device(self.ctx,self.queue,host_matrix)
        self.gpu_result=cl_array.empty(self.queue,(size**2,1,),numpy.float64)
    def combine_kernel(self,arg=""):
        includes="".join(open("includes.cl",'r').readlines())
        defines="".join(open("defines.cl",'r').readlines())
        helpers="".join(open("helpers.cl",'r').readlines())
        arguments="double ix(int i) {double arr[]={"+arg+"}; return arr[i];}"
        program_string=\
            includes+"\n"+\
            defines+"\n"+\
            arguments+"\n"+\
            helpers+"\n"+\
            self.potential+"\n"+\
            self.method
        self.kernel=ElementwiseKernel(self.ctx, "int *x, double start, double end, int size, double *res", \
            "res[i]=get_element(x[i],start,end,size)", preamble=program_string)
    def execute_kernel(self, start, end):
        self.kernel(self.gpu_matrix,start,end,self.size,self.gpu_result)
    def get_results(self):
        return self.gpu_result.get()
