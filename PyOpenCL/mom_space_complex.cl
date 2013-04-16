// The function who's being integrated over.
float2 integrand(float2 r, float2 k, float2 k_prime)
{
    int l=1;
    float j=1.5;
    return c_mul(c_multiplication(r, r, j_l(l, c_mul(k,r)), j_l(l, c_mul(k_prime,r))) , V(r, j));
}
// Calculates one matrix element in berggren-methodology.
float2 get_element_berggren(float2 start, float2 end, float2 weight, float2 xx, float2 yy)
{
    start=c_real(0.0f); // Needed
    end=c_real(7.0f);     // Needed
    float mass=0.019272f;
    float2 sum=to_r(0.0f);
    for (int i=0;i<25;i++)
        sum = sum + c_mul(to_r(get_gauss_legendre_w_25(i)) , (integrand(((end-start)/2.0f*get_gauss_legendre_x_25(i) + (start+end)/2.0f), (xx), (yy))));
    float2 integral=c_mul((end-start),sum)/2.0f;
    // return integral;
    float2 diagonal=c_mul(xx,xx)/(2*mass)*(eps(xx-yy));
    // return diagonal;
    return diagonal+c_mult(2.0f*c_mul(yy,yy),weight,integral/PI);
}
float2 get_element(int x, float2 start, float2 end, int n)
{
    // We got k,k_prim,step
    start=c_real(start);
    end=c_real(10.0f);
    float mass=0.019272f;
    float2 sum=to_r(0.0f);
    float2 step=(end-start)/((float)n);
    float2 xx=c_mul(to_r((float)(x%n)),step);
    float2 yy=c_mul(to_r((float)((x-x%n)/n)),step);
    for (int i=0;i<25;i++)
        sum = sum + c_mul(to_r(get_gauss_legendre_w_25(i)) , (integrand(((end-start)/2.0f*get_gauss_legendre_x_25(i) + (start+end)/2.0f), (xx), (yy))));
    float2 integral=c_mul((end-start),sum)/2.0f;
    // return integral;
    float2 diagonal=c_mul(xx,xx)/(2*mass)*(equal(xx,yy));
    // float2 diagonal=c_mul(xx,xx)/(2*mass)*((x%n)==((x-x%n)/n));
    return diagonal;
    return diagonal+c_mult(2.0f*c_mul(yy,yy),step,integral/PI);
}
// float2 integrand(float2 r, float2 k, float2 k_prime, float2 step)
// {
    // int l=1;
    // float j=1.5;
    // return V(r,j);
    // return c_mul(c_multiplication(r, r, to_r(j_l(l, real(c_mul(k,r)))), to_r(j_l(l, real(c_mul(k_prime,r))))) , V(r, j));
// }

// float integrand(float r, float k, float k_prime, float step)
// {
    // int l=1;
    // float j=1.5;
    // return r*r * j_l(l, k*r) * j_l(l, k_prime*r) * V(r, j);
// }

// float2 get_element(int x, float2 start, float end, int n)
// {
    //return 1.0f;
    // float mass=0.019272f;
    // float2 sum=to_r(0.0f);
    // float2 step=to_r(end-start)/((float)n);
    // float2 xx=((float)(x%n))*step;
    // float2 yy=((float)((x-x%n)/n))*step;
    
    // for (int i=0;i<25;i++)
        // sum = sum + c_mul(to_r(get_gauss_legendre_w_25(i)) , to_r(integrand(real(to_r(end-start)/2.0f*get_gauss_legendre_x_25(i) + to_r(start+end)/2.0f), real(xx), real(yy), real(step))));
        
    // float2 integral=c_mul(to_r(end-start),sum)/2.0f;
    // float2 diagonal=c_mul(xx,xx)/(2*mass)*((x%n)==((x-x%n)/n));
    // return diagonal+c_mult(2.0f*c_mul(yy,yy),step,integral/PI);
// }