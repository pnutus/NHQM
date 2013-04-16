///////////////////
// Complex operations implemented for float2, x-coord for real part and y-coord for imaginary.
//////////////////
float2 c_add(float2 a, float2 b) {return (float2)(a.x+b.x,a.y+b.y);}
float2 c_sub(float2 a, float2 b) {return (float2)(a.x-b.x,a.y-b.y);}
float2 c_mul(float2 a, float2 b) {return (float2)(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);}
float2 c_sin(float2 z) {return (float2)(sin(z.x)*cosh(z.y),cos(z.x)*sinh(z.y));}
float2 c_cos(float2 z) {return (float2)(cos(z.x)*cosh(z.y),sin(z.x)*sinh(z.y));}
float2 c_exp(float2 z) {return (float2)(exp(z.x)*cos(z.y),exp(z.x)*sin(z.y));}
float2 c_mult(float2 a, float2 b, float2 c) {return c_mul(a,c_mul(b,c));}
float2 c_multiplication(float2 a, float2 b, float2 c, float2 d) {return c_mul(c_mul(a,b),c_mul(c,d));}
float2 to_r(float x) {return (float2)(x,0.0f);}
float2 to_i(float y) {return (float2)(0.0f,y);}
float2 to_c(float x, float y) {return (float2)(x,y);}
bool apx_eql(float a, float b, int ulp) {return fabs(a-b) <= ulp*FLT_EPSILON*max(fabs(a), fabs(b));}
bool eps_r(float2 z) {return apx_eql(z.x,0.0f,3);} // Is real part zero?
bool eps_i(float2 z) {return apx_eql(z.y,0.0f,3);} // Is imaginary part zero?
bool eps_c(float2 z) {return eps_r(z)&&eps_i(z);} // Is the complex number zero?
bool eps(float2 z) {return eps_c(z);}
bool equal(float2 a, float2 b) {return apx_eql(a.x,b.x,3)&&apx_eql(a.x,b.x,3);}
float real(float2 z) {return z.x;}
float imag(float2 z) {return z.y;}
float2 c_real(float2 z) {return to_r(z.x);}
float2 c_imag(float2 z) {return to_i(z.y);}
float2 c_conj(float2 z) {return (float2)(z.x,-z.y);}
float abssqr(float2 z) {return z.x*z.x+z.y*z.y;}
float2 epsify(float2 z)
{
    // Is the real part zero? If it is, set it to epsilon.
    // Is the imaginary part zero? Make shure..
    // Useful for divisions of kind sin(x)/x.
    if (eps_r(z))
        z.x=0.0000000001f;
    if (eps_i(z))
        z.y=0.0f;
    return z;
}
float2 c_div(float2 a, float2 b)
{
    float denominator=b.x*b.x+b.y*b.y;
    float2 factor=c_mul(a,c_conj(b));
    return (float2)(factor.x/denominator,factor.y/denominator);
}