// The Woods-Saxon potential.
float2 V(float2 r, float j)
{
    float V0=-47.0f, Vso=-7.5f, r0=2.0f, d=0.65f, l=1.0f;
    float2 f=c_div(to_r(1.0f),(to_r(1.0f)+c_exp((r-to_r(r0))/d)));
    float spin_orbit=0.5f*(j*(j+1.0f)-l*(l+1.0f)-0.75f);
    return c_mul(f,  to_r(V0)-c_mul(to_r(4.0f*Vso*spin_orbit)  ,  c_div((f-to_r(1.0f)) , c_mul(to_r(d),r))));
}
