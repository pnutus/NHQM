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

kernel=ElementwiseKernel(ctx,"float *num, int items, float2 *eigvecs, float2* weights, float2* points, float2 *res", "res[i]=func(num[i], eigvecs)",\
                            preamble="""
                            float func(float a, __global float *b)
                            {
                            int i;
                            float2 inner_sum = 0;
                            for(i=0; i < items; i++)
                            {
                                inner_sum += wights[i] , eigvecs[i] , ********
                            }
                            
                            }
                            
                            float2 V(float2 r)
                            {
                            	float2 beta = to_r(1.0f);
                            	float2 res = c_exp(EULER , c_mult(-beta, r));
                            	return res;	
                            }

                            float2 integrand(float2 r, float2 k, float2 k_prim)
                            {
                                int l=1;
                                float j=1.5;
                                return c_mul(c_multiplication(r, r, j_l(l, c_mul(k,r)),
                                j_l(l, c_mul(k_prime,r))) , V(r, j));
                            }
                            
                            float2 v_sep(float2 k, float2 k_prim, int items)
                            {
                                int i;
                                float2 res = 0;
                                for (i=0; i< items; i++)
                                {
                                    res += integrand(weights[i], k, k_prim)
                                }
                                return c_multiplication(to_r(0.5f), to_r(PI), c_pow(k_prim, 2), res)     
                            }
                            
                            
                            """)
kernel(gpu_arr,gpu_arr,gpu_res)
print gpu_res.get()


# Does not work currently...
# It should return 0,2,4,6, etc, it doesn't