#include <cstdio>
#include <iostream>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/ZZ_pE.h>

using namespace NTL;
//s - Vector of L polynomials From Rq (same as Rp but with q as the p)
//Rp set of some polynomial with coefficents less than P. P is a the modulus (is infinitly large)
// b' Vector of polynomials in RP
//s % p - 
//R2 set of polynomials coefficents are less than 2 for example, 1 + x^2 + x^100     is in R_2

//Eq = 12
//Ep = 10
//Et = 3

#define EQ 12
#define EP 10
#define ET 3

int main(void){
    /* Line 1 */
    
    // b' s and c come from other steps 
    int l = 10; //TODO - figure out what this should be.
    constexpr int degree = 256; //TODO - Change to 256 Later
    ZZ p(1024);
    ZZ_p::init(p);

    //ZZ constPow(128); // 2 ^ 7
    ZZ_pX constPow(128);
    std::cout<<constPow;

    ZZ pow1; //=
    ZZ pow2; //=
    ZZ pow3; //=
    power(pow1, 2,EP-1);
    power(pow2,2,EP - ET - 1);
    power(pow3,2,EQ - EP - 1);
    
    ZZ T(8);
    
    const ZZ h2Val = pow1 - pow2 + pow3;
    ZZ_p h2Val_p;// = conv(h2Val);
    conv(h2Val_p, h2Val);
    ZZ_pX h2;

    for(size_t i = 0 ; i<degree;i++){
        SetCoeff(h2,i,h2Val_p);
    }

    ZZ q;
    q = ZZ(4096);
    

    ZZ_pX polyMod;
    polyMod.SetLength(degree + 1);
    SetCoeff(polyMod, 0,1);
    SetCoeff(polyMod, degree,1);

    ZZ_pX v;
    v.SetLength(8);

    Vec<ZZ_pX> s;
    Vec<ZZ_pX> b;
    
    s.SetLength(l);
    b.SetLength(l);

    for(size_t i = 0 ; i < s.length() ; i++){
        random(s[i], degree);
        random(b[i], degree);
    }

    //b'T (s % p)
    for(size_t i = 0; i < l ; i++){
        v+= s[i] * b[i];
    }

    //  s % p
    v = v % polyMod;

    ZZ_p::init(ZZ(8)); // 8 is T
    ZZ_pX cm; // Get this from encryption (cypher text)
    random(cm,degree);
    std::cout << '\n';


    ZZ_p::init(p);
    ZZ_pX temp1;

    temp1 = v - constPow * cm + h2;

    ZZ_pX temp2;

    for (size_t i = 0 ; i < degree; i++){
       // temp2[i] = temp1[i] % p;
       SetCoeff(temp2,i,ZZ_p(temp1[i]));
    }

    std::cout << temp2 <<'\n';
    for(size_t i = 0 ; i < degree ; i++){
        unsigned int temp3;
        conv(temp3, temp2[i]);
        temp3 >>= (EP-1);
        temp2[i] = temp3;
        //temp2[i] /= 512; // >> EP-1
    }
    std::cout << temp2 <<'\n';
    
    return 0; 
    
}