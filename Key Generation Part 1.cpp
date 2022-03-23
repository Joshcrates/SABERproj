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

#include "crypto++/cryptlib.h"
#include "crypto++/shake.h"
#include "crypto++/filters.h"
using CryptoPP::SHAKE128;


using namespace std;

int main(int argc, char** argv) {

    //Step 1 - seed_A ← U({0,1}^256)        Generates a 256 byte seed consisting of 1s and 0s (Uniform).

    //setup to generate uniform distribution of 1s and 0s
    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<> distrib(0, 1);

    //generate seed_A
    vector<int> seed_A;
    //loop to generate 256 bytes
    for(int i = 0; i < 256; i++) {
        //loop to generate 8 bits for each byte
        for(int j = 0; j < 8; j++) {
            seed_A.push_back(distrib(gen));
        }
    }


    //Step 2 - A = gen(seed_A) ∈ R_q^l×l    Generates a pseudorandom matrix using seed A and SHAKE-128.

    //Using SHAKE128 from crypto++, this only generates a SHAKE128 hash, and not a psuedorandom matrix.
    std::string digest;

    SHAKE128 hash;
    //something wrong with this line below:
    hash.Update((const CryptoPP::byte*)seed_A.data(), seed_A.size());
    digest.resize(hash.DigestSize());
    hash.Final((CryptoPP::byte*)&digest[0]);

    CryptoPP::StringSource(digest, true);   //double check this part for proper arguments (encoder is currently default value)

    //using the shake128 string, a matrix sized lxl is populated and set to variable A


    //Step 3 - r = U({0,1}^256)     Generates a 256 byte key consisting of 1s and 0s (Uniform distribution).
    //The process is the same as step 1.

    vector<int> r;
    //loop to generate 256 bytes
    for(int i = 0; i < 256; i++) {
        //loop to generate 8 bits for each byte
        for(int j = 0; j < 8; j++) {
            r.push_back(distrib(gen));
        }
    }


    //test seed generation
    cout << "[";
    for(int i = 0; i < r.size(); i++) {
        cout << r[i] << ", ";
    }
    cout << "]" << endl;


    //test shake128
    std::cout << "Digest: ";
    for (unsigned char c : digest){
            printf("%x",c);
            }
    std::cout << std::endl;


    return 0;
}


//g++ -I/usr/local/include myprog1.cpp -o myprog1.out -L/usr/local/lib -lcryptopp
//./myprog1.out
