//Key Generation Part 1 - James Goedmakers

#include <cstdlib>
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <sstream>
#include <numeric>
#include <cmath>
#include <random>

#include "cryptopp-master/shake.h"
using CryptoPP::SHAKE128;


using namespace std;

int main(int argc, char** argv) {

    //Step 1 - seed_A ← U({0,1}^256)        Generates a 256 bit seed consisting of 1s and 0s (Uniform).

    //setup to generate uniform distribution of 1s and 0s
    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<> distrib(0, 1);

    //generate seed_A
    vector<int> seed_A;
    for(int i = 0; i < 256; i++) {
        seed_A.push_back(distrib(gen));
    }
    

    //Step 2 - A = gen(seed_A) ∈ R_q^l×l    Generates a pseudorandom matrix using seed A and SHAKE-128.
    //Using SHAKE128 from crypto++, this only generates a SHAKE128 hash, and not a psuedorandom matrix.
    //A = SHAKE128(seed_A)
    
    // TODO: This causes a linker error on MSVC (Need to test on GCC and clang)
    //SHAKE128 hash;
    //hash.Update((const CryptoPP::byte*)seed_A.data(), seed_A.size());

    //Step 3 - r = U({0,1}^256)     Generates a 256 bit key consisting of 1s and 0s (Uniform distribution).
    //The process is the same as step 1.
    vector<int> r;
    for(int i = 0; i < 256; i++) {
        r.push_back(distrib(gen));
    }

    //test 
    cout << "[";
    for(int i = 0; i < r.size(); i++) {
        cout << r[i] << ", ";
    }
    cout << "]" << endl;

    return 0;
}