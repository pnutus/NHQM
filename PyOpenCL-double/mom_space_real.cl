double integrand(double r, double k, double k_prime, double step)
{
    int l=1;
    double j=1.5;
    return r*r * j_l(l, k*r) * j_l(l, k_prime*r) * V(r, j);
}

double get_element(int x, double start, double end, int n)
{
    double mass=0.019272;
    double* w,z;
    double sum=0.0;
    double step=(end-start)/((double)n);
    double xx=((double)(x%n))*step;
    double yy=((double)((x-x%n)/n))*step;
    for (int i=0;i<25;i++)
    {
        sum += (get_gauss_legendre_w_25(i) * integrand( ((end-start)/2.0)*get_gauss_legendre_x_25(i) + (start+end)/2.0, xx, yy, step));
    }
    double integral=(end-start)*sum/2.0;
    double diagonal=(xx*xx)/(2*mass)*((x%n)==((x-x%n)/n));
    return diagonal+2*(yy*yy)*step*integral/PI;
}
