double j_l(int l, double r)
{
    switch (l)
    {
        case 0:
            if (fabs(r)<0.0001) return 1.0;
            else return sin(r)/r;
        case 1:
            if (fabs(r)<0.0001) return 0.0;
            else return (sin(r)/(r*r)-cos(r)/r);
        case 2:
            return ((3/(r*r)-1)*sin(r)/r-3*cos(r)/(r*r));
        case 3:
            return (-(15/(r*r*r)+6/r)*sin(r)/r-(15/(r*r)-1)*cos(r)/r);
    }
}
