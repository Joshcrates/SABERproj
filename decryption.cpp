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

int main0(void){
    ZZ p(7);
    ZZ_p::init(p);

    long degree = 16;

    ZZ_pX testPoly;
    testPoly.SetLength(16);

    SetCoeff(testPoly, 3,25); // 25 % 7 = 4

    Vec<ZZ_pX> vector;
    vector.SetLength(5);

    ZZ_pX temp;

    for(size_t i = 0 ; i < vector.length() ; i++){
        random(vector[i], degree);
        std::cout << vector[i] << '\n';
    }

return 0; 
   // std::cout <<testPoly;
}

int main(void){
    /* Line 1 */
    
    // b' s and c come from other steps 
    int l = 10; //TODO - figure out what this should be.
    constexpr int degree = 256; //TODO - Change to 256 Later
    ZZ p(1024);
    ZZ_p::init(p);

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

    std::cout <<"polyMod    " << polyMod << '\n';
    std::cout << v <<'\n';
    
    return 0; 
    
}