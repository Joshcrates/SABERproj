//Encryption Part 1 - Joshua Walsworth

#include <cstdlib>
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <sstream>
#include <numeric>
#include <cmath>
#include <random>

#include "shake.h"
using CryptoPP::SHAKE128;


using namespace std;

int main(int argc, char** argv) {


    //setup to generate uniform distribution of 1s and 0s
    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<> distrib(0, 1);

    //This should be passed as an argument but our project is not organized in such a way for that to be possible for now.
    //generate seed_A
    vector<int> seed_A;
    for (int i = 0; i < 256; i++) {
        seed_A.push_back(distrib(gen));
    }

    //Note: Since this line is identical to line 2 of keygen, I am using the code James wrote.
    //Step 1 - A = gen(seed_A) ∈ R_q^l×l    Generates a pseudorandom matrix using seed A and SHAKE-128.
    //Using SHAKE128 from crypto++, this only generates a SHAKE128 hash, and not a psuedorandom matrix.
    //A = SHAKE128(seed_A)

    /* Commenting out this line because it does not compile and thus needs further work
    SHAKE128 hash;
    //hash.Update((const byte*)seed_A.data(), seed_A.size());
    */

    //Step 2 - if r is not specified then:      Checks if r has been passed as an argument.
    //Due to the incomplete nature of our project at this moment, our code is not organized enough 
    //to have separate functaions all on one file. For that reason, 
    //I am going to assume for now that r is not passed as an argument.

    //Note: Again, this line is the same same as the line 3 James wrote for keygen
    //Step 3 - r = U({0,1}^256)     Generates a 256 bit key consisting of 1s and 0s (Uniform distribution).
    //The process is the same as step 1.
    vector<int> r;
    for (int i = 0; i < 256; i++) {
        r.push_back(distrib(gen));
    }

    return 0;
}

void alg15() {
    int l = 2;
    int eq = 12;
    int n = 256;
    int length = l * l * n * eq / 8;
    cout << length << endl;
    CryptoPP::byte buf[1536];//1024*12/8 = 1536*8 = 12288 total bits
    for(int i=0;i<1536;i++)
        buf[i] = '0';
    cout << buf[0] << buf[length - 1];

    //goal: understand byte data type so I can sort them into bits
    //goal: learn SHAKE-128 syntax to implement alg15 line 2
    bool bufs[1024][12];//1024*12 = 12288 total bits
    for (int i = 0; i < 1024; i++) {
        bufs[i]=
    }

}

