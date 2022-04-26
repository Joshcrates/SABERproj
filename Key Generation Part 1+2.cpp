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
#include <bitset>

#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/ZZ_pE.h>

#include "cryptopp860/shake.h"
#include "cryptopp860/cryptlib.h"
#include "cryptopp860/filters.h"
using CryptoPP::SHAKE128;

const int l = 2;
const int n = 256;
const int eq = 12;  //check later if 13 or 12

using namespace std;
using namespace NTL;

int main(int argc, char** argv) {

    //Step 1 - seed_A ← U({0,1}^256)
    //Generates a 256 byte seed consisting of 1s and 0s (Uniform).

    //setup to generate uniform distribution of 1s and 0s
    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<> distrib(0, 1);

    //generate seed_A
    vector<int> seed_A;
    //loop to generate 256 bytes, DOUBLE CHECK THIS- SABER_SEEDBYTES SAYS 32
    for (int i = 0; i < 256; i++) {
        //loop to generate 8 bits for each byte
        for (int j = 0; j < 8; j++) {
            seed_A.push_back(distrib(gen));
        }
    }


    //Step 2 - A = gen(seed_A) ∈ R_q^l×l    
    //Generates a pseudorandom matrix using seed A and SHAKE-128.

    //Using SHAKE128 from crypto++, this only generates a SHAKE128 hash, and not a psuedorandom matrix.
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

    //Step 3 - r = U({0,1}^256)     
    //Generates a 256 byte key consisting of 1s and 0s (Uniform distribution).

    //The process is the same as step 1.
    vector<int> r;
    //loop to generate 256 bytes
    for (int i = 0; i < 256; i++) {
        //loop to generate 8 bits for each byte
        for (int j = 0; j < 8; j++) {
            r.push_back(distrib(gen));
        }
    }

    //test seed generation
    /*cout << "[";
    for(int i = 0; i < r.size(); i++) {
        cout << r[i] << ", ";
    }
    cout << "]" << endl;*/





    /*
    ZZ_pXMatrix poly;

    ZZ size = ZZ(1<<13);
    ZZ_p::init(size);

    vector<ZZ_pX> s;

    int i, j;
    zz_pX z;
    */

    //step 4 (NTL was acting weird so I just used an array to represent the polynomial)
    default_random_engine gen2(0);
    binomial_distribution<int> binomial(8, 0.5);

    int i, j, s[16][16][3];
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            for (k = 0; k < 3; k++)
                s[i][j][k] = ((binomial(gen2) - 4) % (1 << 13));
        }
    }
    //testing step 4
    /*
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            for (k = 0; k < 3; k++) {
                cout << s[i][j][k] << ' ';
            }
        }
        cout << endl;
    }
    */
    //step 5
    // I also had to generate A randomly because step 2 is not fully implemented
    int A[16][16][64], newA[16][16][64] = { 0 }, temp[16][16][64] = { 0 };
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            for (k = 0; k < 3; k++) {
                A[i][j][k] = (binomial(gen2) % (1 << 13));
                newA[i][j][k] = A[i][j][k];
            }
        }
    }

    int h, i2, j2, k1, k2;
    //multiply A by itself 2^4 times and then by s (should probably get NTL working first)
    for (h = 0; h < 15; h++) {
        //multiply matrices
        for (i = 0; i < 16; i++)
            for (j = 0; j < 16; j++)
                for (i2 = 0; i2 < 16; i2++)
                    for (j2 = 0; j2 < 16; j2++)
                        for (k1 = 0; k1 < 32; k1++)
                            for (k2 = 0; k2 < 3; k2++) {
                                temp[i][j][k1 + k2] += newA[i2][j][k1] * A[i][j2][k2];
                                temp[i][j][k1 + k2] %= 1 << 13;
                            }

        //set newA to temp and clear temp
        for (i = 0; i < 16; i++)
            for (j = 0; j < 16; j++)
                for (k = 0; k < (3 + 2 * h); k++) {
                    newA[i][j][k] = temp[i][j][k];
                    temp[i][j][k] = 0;
                }
    }
    //multiply by s
    for (i = 0; i < 16; i++)
        for (j = 0; j < 16; j++)
            for (i2 = 0; i2 < 16; i2++)
                for (j2 = 0; j2 < 16; j2++)
                    for (k1 = 0; k1 < 32; k1++)
                        for (k2 = 0; k2 < 3; k2++) {
                            temp[i][j][k1 + k2] += newA[i2][j][k1] * s[i][j2][k2];
                            temp[i][j][k1 + k2] %= 1 << 13;
                        }
    for (i = 0; i < 16; i++)
        for (j = 0; j < 16; j++)
            for (k = 0; k < 64; k++) {
                newA[i][j][k] = temp[i][j][k];
                temp[i][j][k] = 0;
            }

    //add h to the result and shift to the right 3 bits
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            for (k = 0; k < 27; k++) {
                newA[i][j][k] += 4;
                newA[i][j][k] %= (1 << 13);
                newA[i][j][k] >>= 3;
            }
        }
    }

    //test results
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            cout << i << ',' << j << ": ";
            for (k = 34; k > 0; k--) {
                if (newA[i][j][k])
                    cout << newA[i][j][k] << "x^" << k - 1 << ' ';
            }
            cout << endl;
        }
        cout << endl;
    }

    return 0;
}