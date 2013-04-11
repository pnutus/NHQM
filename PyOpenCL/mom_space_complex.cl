float integrand(float2 r, float2 k, float2 k_prime, float2 step)
{
    int l=1;
    float j=1.5;
    return r*r * j_l(l, k*r) * j_l(l, k_prime*r) * V(r, j);
}

float get_element(int x, float2 k, float2 k_prim, float2 w, float2 start, float2 end, int n)
{
    float mass=0.019272;
    float* w,z;
    float sum=0.0;
    float step=(end-start)/((float)n);
    float xx=((float)(x%n))*step;
    float yy=((float)((x-x%n)/n))*step;
    for (int i=0;i<25;i++)
    {
        sum += (get_gauss_legendre_w_25(i) * integrand( ((end-start)/2.0)*get_gauss_legendre_x_25(i) + (start+end)/2.0, xx, yy, step));
    }
    float integral=(end-start)*sum/2.0;
    float diagonal=(xx*xx)/(2*mass)*((x%n)==((x-x%n)/n));
    return diagonal+2*(yy*yy)*step*integral/PI;
}