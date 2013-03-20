float integrand(float r) {return V(r,1);}
float get_element(int x, float start, float end, int n)
{
    //return get_gauss_legendre_w_25(3);
    float sum=0.0;
    for (int i=0;i<25;i++)
        sum += (get_gauss_legendre_w_25(i) * integrand( ((end-start)/2.0)*get_gauss_legendre_x_25(i) + (start+end)/2.0));
    return sum*((end-start)/((float)n));
}