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
    def set_method(self,method="mom_space_real.cl"):
        self.method="".join(open(method,'r').readlines())
    def allocate_space(self,type,contour):
        self.size=len(contour[1])
        (k,w)=(contour)
        host_matrix=(numpy.array([i for i in range(self.size**2)])).astype(numpy.int32)
        host_k=(numpy.array([k[i//self.size] for i in range(self.size**2)])).astype(type)
        host_k_prim=(numpy.array([k[i%self.size] for i in range(self.size**2)])).astype(type)
        host_w=(numpy.array([w[i%self.size] for i in range(self.size**2)])).astype(type)
        self.gpu_k=cl_array.to_device(self.ctx,self.queue,host_k)
        self.gpu_k_prim=cl_array.to_device(self.ctx,self.queue,host_k_prim)
        self.gpu_w=cl_array.to_device(self.ctx,self.queue,host_w)
        self.gpu_matrix=cl_array.to_device(self.ctx,self.queue,host_matrix)
        self.gpu_result=cl_array.empty(self.queue,(self.size**2,1,),type)
    def combine_kernel(self,arg=""):
        includes="".join(open("includes.cl",'r').readlines())
        defines="".join(open("defines.cl",'r').readlines())
        helpers="".join(open("helpers_complex.cl",'r').readlines())
        arguments="float ix(int i) {float arr[]={"+arg+"}; return arr[i];}"
        program_string=\
            includes+"\n"+\
            defines+"\n"+\
            arguments+"\n"+\
            helpers+"\n"+\
            self.potential+"\n"+\
            self.method
        self.kernel=ElementwiseKernel(self.ctx, "int *x,c_float *k,c_float *k_prim,c_float *w, float start, float end, int size, c_float *res", \
            "res[i]=get_element(x[i],k[i],k_prim[i],w[i],start,end,size)", preamble=program_string)
    def execute_kernel(self, start, end):
        self.kernel(self.gpu_matrix,self.gpu_k,self.gpu_k_prim,self.gpu_w,start,end,self.size,self.gpu_result)
    def get_results(self):
        return self.gpu_result.get()