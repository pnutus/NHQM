float2 get_element(float2 start, float2 end, float2 weight, float2 xx, float2 yy)
{

    // start=c_real(0.0f); // Needed
    // end=c_real(7.0f); // Needed
    // float mass=0.019272f;
    // float2 sum=to_r(0.0f);
    // for (int i=0;i<25;i++)
    //     sum = sum + c_mul(to_r(get_gauss_legendre_w_25(i)) , (integrand(((end-start)/2.0f*get_gauss_legendre_x_25(i) + (start+end)/2.0f), (xx), (yy))));
    // float2 integral=c_mul((end-start),sum)/2.0f;
    // // return integral;
    // 	
    // 
    // float2 diagonal=c_div(c_mult(xx,xx,(eps(xx-yy))),to_r(2*mass));
    // // return diagonal;
    // return diagonal+c_mult(2.0f*c_mul(yy,yy),weight,c_div(integral , to_r(PI) ) );	
}

	
float2 V(float2 r)
{
	float2 beta = to_r(1.0f);
	float2 res = pow(EULER , c_mult(-beta, r));
	// Det finns en funktion c_exp som gör detta.
	// Har du inte skrivit pow själv fixar den tyvärr inte komplexa tal.
	return res;	
}

float2 integrand(float2 r, float2 k, float2 k_prime)
{
    int l=1;
    float j=1.5;
    return c_mul(c_multiplication(r, r, j_l(l, c_mul(k,r)), j_l(l, c_mul(k_prime,r))) , V(r, j));
}


// 
// def potential(r, l, j):
//     return sp.exp(- beta * r**2)

// # def mom.integrand(r, k, k_prim, V, l, j):
// #     return r**2 * j_l(l, k*r) * j_l(l, k_prim*r) * V(r, l, j)
