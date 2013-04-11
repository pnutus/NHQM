double V(double r, double j)
{
    double V0=-70.0, Vso=-7.5, r0=2.0, d=0.65, l=1.0;
    double f=1/(1+exp((r-r0)/d));
    double spin_orbit=0.5*(j*(j+1)-l*(l+1)-0.75);
    return f*(V0-4*Vso*spin_orbit*(f-1)/(d*r));
}
