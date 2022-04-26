//Encryption Part 1 - Joshua Walsworth

#include <cstdlib>
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <sstream>
#include <numeric>
#include <cmath>
#include <bitset>
#include <random>

#include <cryptopp/shake.h>
#include <cryptopp/cryptlib.h>
#include <cryptopp/filters.h>

using CryptoPP::SHAKE128;


const int l = 2;
const int n = 256;
const int eq = 12;

using namespace std;

void alg15();

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
    //Step 1 -> A = gen(seed_A) ∈ R_q^l×l    Generates a pseudorandom matrix using seed A and SHAKE-128.
    //Using SHAKE128 from crypto++, this only generates a SHAKE128 hash, and not a psuedorandom matrix.
    //A = SHAKE128(seed_A)
    //Generates a pseudorandom matrix using seed A and SHAKE-128.

    std::string digest;

    SHAKE128 hash;
    hash.Update((const CryptoPP::byte*)seed_A.data(), seed_A.size());
    digest.resize(hash.DigestSize());
    hash.Final((CryptoPP::byte*)&digest[0]);

    CryptoPP::StringSource(digest, true);   //double check this part for proper arguments (encoder is currently default value)

    //using the shake128 string, a matrix sized lxl is populated and set to variable A:

    //Matrix generation:
    //Instantiate byte string object buf of length l^2 × n × ϵq/8
    vector<char> buf(l * l * n * eq / 8, '0');

    //generate seed_B
    vector<int> seed_B;
    //loop to generate 256 bytes
    for (int i = 0; i < l * l * n * eq / 8; i++) {
        //loop to generate 8 bits for each byte
        for (int j = 0; j < 8; j++) {
            seed_B.push_back(distrib(gen));
        }
    }

    //SHAKE-128(buf, l^2 × n × ϵq/8, seedAAA, SABER SEEDBYTES)
    hash.Update((const CryptoPP::byte*)seed_B.data(), seed_B.size());
    //buf.resize(hash.DigestSize());
    hash.Final((CryptoPP::byte*)&buf[0]);

    CryptoPP::StringSource(&buf[0], true);
    string buf_as_a_string(buf.begin(), buf.end());
    cout << "buf_as_a_string size:" << buf_as_a_string.size() << endl;
    cout << "seed_B size:" << seed_B.size() << endl;

    //Split buf into l^2 × n equal bit strings of bit length ϵq and obtain
    //(bufl2n−1 ‖ . . . ‖ buf0) = buf
    vector<bitset<eq>> buf_splitted(l * l * n, '0');
    for (int i = 0; i < l * l * n; i++) {
        buf_splitted[i] = bitset<eq>(buf_as_a_string.c_str()[i]);
    }

    cout << buf_splitted[0] << endl;
    cout << buf_splitted[l * l * n - 1] << endl;  //check the bit/byte lengths, why are these all zeros?
    //exit(0);
    //the matrix is 2 dimensional, so nested for loops are used
    int matrix_A[l][l][n];
    int k = 0;
    for (int i_1 = 0; i_1 < l; i_1++) {
        for (int i_2 = 0; i_2 < l; i_2++) {
            for (int j = 0; j < n; j++) {
                //printf("x");
                //matrix_A[i_1][i_2][j] = 1;
                matrix_A[i_1][i_2][j] = buf_splitted[k].to_ulong();
                k++;
            }
        }
    }

    for (int i_1 = 0; i_1 < l; i_1++) {
        for (int i_2 = 0; i_2 < l; i_2++) {
            cout << matrix_A[i_1][i_2] << " ";
        }
        cout << endl;
    }

    //Step 2 - if r is not specified then:      Checks if r has been passed as an argument.
    //Due to the incomplete nature of our project at this moment, our code is not organized enough 
    //to have separate functions all on one file. For that reason, 
    //I am going to assume for now that r is not passed as an argument.

    //Note: Again, this line is the same as the line 3 James wrote for keygen
    //Step 3 - r = U({0,1}^256)     Generates a 256 bit key consisting of 1s and 0s (Uniform distribution).
    //The process is the same as step 1.
    vector<int> r;
    //loop to generate 256 bytes
    for (int i = 0; i < 256; i++) {
        //loop to generate 8 bits for each byte
        for (int j = 0; j < 8; j++) {
            r.push_back(distrib(gen));
        }
    }

    alg15();

    return 0;
}

void alg15() {
    int l = 2;
    int eq = 12;
    int n = 256;
    int length = l * l * n * eq / 8;
    cout << length << endl;
    CryptoPP::byte buf[1536];//1024*12/8 = 1536*8 = 12288 total bits
    for (int i = 0; i < 1536; i++)
        buf[i] = '0';
    cout << buf[0] << buf[length - 1];

    //goal: understand byte data type so I can sort them into bits
    //goal: learn SHAKE-128 syntax to implement alg15 line 2
    bool bufs[1024][12];//1024*12 = 12288 total bits
    for (int i = 0; i < 1024; i++) {
        //      bufs[i]=
    }

}
