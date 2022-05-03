//SABER Key Generation
//Part 1 - James Goedmakers
//Part 2 - Robert Blanchette

//SABER Encryption
//Joshua Walsworth

//SABER Decryption
//Louis Traficante

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

//NTL Library used for polynomials
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/ZZ_pE.h>

//Crypto++ library used for SHAKE128
#include <crypto++/cryptlib.h>
#include <crypto++/shake.h>
#include <crypto++/filters.h>
#include <crypto++/hex.h>
#include <crypto++/files.h>

using CryptoPP::SHAKE128;

//Default SABER constants
const int l = 3;
const int n = 256;
const int eq = 13; 
const int EP = 10;
const int ET = 4;
const int p = pow(2,10);
const int q = pow(2,13);
const int T = pow(2,ET); //table 2 page 12

#define SABER_SEEDBYTES 32;

using namespace std;

int main(int argc, char** argv) {
    using namespace CryptoPP;
    using namespace NTL;

    //******************************KEYGEN******************************
    //******************************STEP 1******************************
    //Step 1 - seed_A ← U({0,1}^256)
    //Generates a 256 byte seed consisting of 1s and 0s (Uniform).
    //setup to generate uniform distribution of 1s and 0s
    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<> distrib(0, 1); //uniform distribution as specified by the algorithm

    ZZ mod_p = ZZ(7);
    ZZ_p::init(mod_p);

    //generate seed_A
    vector<int> seed_A;
    //loop to generate 256 bytes, DOUBLE CHECK THIS- SABER_SEEDBYTES SAYS 32
    for(int i = 0; i < 256; i++) {
        //loop to generate 8 bits for each byte
        for(int j = 0; j < 8; j++) {
            seed_A.push_back(distrib(gen));
        }
    }

    //******************************STEP 2******************************
    //Step 2 - A = gen(seed_A) ∈ R_q^l×l    
    //Generates a pseudorandom matrix using seed A and SHAKE-128.
    //using the shake128 string, a matrix sized lxl is populated and set to variable A:

    //Instantiate byte string object buf of length l^2 × n × ϵq/8
    //vector<char> buf(l*l*n*eq/8,'0');
    std::string buf;

    //Using SHAKE128 from crypto++, this only generates a SHAKE128 hash, and not a psuedorandom matrix.
    SHAKE128 hash;
    //SHAKE-128(buf, l^2 × n × ϵq/8, seedAAA, SABER SEEDBYTES)
    hash.Update((const byte*)seed_A.data(), seed_A.size());
    //buf.resize(hash.DigestSize());
    buf.resize(32);
    hash.Final((CryptoPP::byte*)&buf[0]);
    //CryptoPP::StringSource(&buf[0], true);
    //string buf_as_a_string(buf.begin(),buf.end());
    //CryptoPP::HexEncoder encoder;
    CryptoPP::HexEncoder encoder(new FileSink(std::cout));
    //StringSource(buf, true, new Redirector(encoder));

    //Split buf into l^2 × n equal bit strings of bit length ϵq and obtain
    //(bufl2n−1 ‖ . . . ‖ buf0) = buf
    //vector<bitset<eq>> buf_splitted(l*l*n,'0');
    vector<bitset<eq>> buf_splitted;
    for(int i = 0; i < l*l*n; i++) {
        buf_splitted.push_back(bitset<eq>(buf[i]));
    }

    //Matrix generation:
    //the matrix is 2 dimensional, so nested for loops are used
    //matrix_A is filled with buf splits
    int matrix_A[l][l][n];
    int k = 0;
    for(int i_1 = 0; i_1 < l; i_1++) {
        for(int i_2 = 0; i_2 < l; i_2++) {
            for(int j = 0; j < n; j++) {
                matrix_A[i_1][i_2][j] = buf_splitted[k].to_ulong();
                k++;
            }
        }
    }    

    Mat<ZZ_pX> matrix_A_NTL;
    matrix_A_NTL.SetDims(l,l);
    k = 0;
    for(int i_1 = 0; i_1 < l; i_1++) {
        for(int i_2 = 0; i_2 < l; i_2++) {
            for(int j = 0; j < n; j++) {
                //matrix_A[i_1][i_2][j] = buf_splitted[k].to_ulong();
                matrix_A_NTL[i_1][i_2].SetLength(l);
                SetCoeff(matrix_A_NTL[i_1][i_2],j,ZZ_p(buf_splitted[k].to_ulong()));
                k++;
            }
        }
    }    

    //cout << "Testing matrix A NTL: " << endl;
    //cout << matrix_A_NTL[1][1] << endl;
    //exit(0);

    //******************************STEP 3******************************
    //Step 3 - r = U({0,1}^256)     
    //Generates a 256 byte key consisting of 1s and 0s (Uniform distribution).
    //The process is the same as step 1.
    vector<int> r;
    //loop to generate 256 bytes
    for(int i = 0; i < 256; i++) {
        //loop to generate 8 bits for each byte
        for(int j = 0; j < 8; j++) {
            r.push_back(distrib(gen));
        }
    }

    //******************************STEP 4******************************

    //s = βμ(Rql×1 ; r)

    ZZ p2 = ZZ(7);
    ZZ_p::init(p2);
    NTL::ZZ_pXMatrix poly;
    NTL::ZZ size = NTL::ZZ(1<<13);
    NTL::ZZ_p::init(size);
    Vec<NTL::ZZ_pX> s2;
    s2.SetLength(l);
    NTL::zz_pX z;
    

    //step 4 (NTL was acting weird so I just used an array to represent the polynomial)
    default_random_engine gen2(0);
    binomial_distribution<int> binomial(8, 0.5);

    int i, j, s[16][16][l];
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            for (k = 0; k < l; k++)
                s[i][j][k] = ((binomial(gen2) - 4) % (1 << eq));
        }
    }

    k = 0;
    for (i = 0; i < l; i++) {
        for (j = 0; j < n; j++) {
            //we need to implement s[i][j] = HammingWeight(bufk) − HammingWeight(bufk+1) mod q
            ZZ_p temp = ZZ_p(buf_splitted[k].count()-buf_splitted[k+1].count()); //check if count returns hamming weight of bitset
            SetCoeff(s2[i],j,temp); //s2 is the final secret key
            k = k + 2;
        }
    }

    //cout << "Testing s2[l-1]" << endl;
    //cout << s2[l-1] << endl;
    //exit(0);

    
    //******************************STEP 5******************************
    //b = ((A^T s + h) mod q) >> (eq − ep) ∈ Rpl×1

    ZZ_p::init(ZZ(q));

    //Spec page 8: vector h is made up of h1
    const ZZ h1Val = ZZ(eq - EP - 1);
    ZZ_p h1Val_p;// = conv(h2Val);
    conv(h1Val_p, h1Val);
    ZZ_pX h1;

    for(size_t i = 0 ; i<n;i++){ //n = degree
        SetCoeff(h1,i,h1Val_p);
    }

    Vec<ZZ_pX>h_vec;
    h_vec.SetLength(l);

    for(size_t i = 0 ; i<l;i++){
        h_vec[i] = h1;
    }

    Vec<ZZ_pX> b_vec;
    b_vec.SetLength(l);

    for (i = 0; i < l; i++) {
        for (j = 0; j < l; j++) {
            b_vec[i] = matrix_A_NTL[j][i]*s2[i] + h_vec[i];
        }
    }

    for (i = 0; i < l; i++) {
        for (j = 0; j < l; j++) {
            //b_vec[i] >> (eq - EP);
            unsigned int temp3;
            conv(temp3, b_vec[i][j]);
            temp3 >>= (eq - EP);
            b_vec[i][j] = temp3;
        }
    }

    //step 5
    // I also had to generate A randomly because step 2 is not fully implemented
    int A[16][16][64], newA[16][16][64] = { 0 }, temp[16][16][64] = { 0 };
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            for (k = 0; k < l; k++) {
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
                            for (k2 = 0; k2 < l; k2++) {
                                temp[i][j][k1 + k2] += newA[i2][j][k1] * A[i][j2][k2];
                                temp[i][j][k1 + k2] %= 1 << 13;
                            }

        //set newA to temp and clear temp
        for (i = 0; i < 16; i++)
            for (j = 0; j < 16; j++)
                for (k = 0; k < (l + 2 * h); k++) {
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
                        for (k2 = 0; k2 < l; k2++) {
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
                newA[i][j][k] >>= l;
            }
        }
    }

    //Vec<NTL::ZZ_pX> b;
    //b.SetLength(l);
    
    //******************************STEP 6******************************
    //return (pk := (seedA, b), s) (seed_A, b_vec, s2)

    //****************************ENCRYPTION****************************
    
    /**
    int length = l * l * n * eq / 8;
    //cout << length << endl;
    CryptoPP::byte buf3[1536];//1024*12/8 = 1536*8 = 12288 total bits
    for (int i = 0; i < 1536; i++) {
        buf3[i] = '0';
    }
    //cout << buf3[0] << buf3[length - 1];

    //goal: understand byte data type so I can sort them into bits
    //goal: learn SHAKE-128 syntax to implement alg15 line 2
    bool bufs[1024][12];//1024*12 = 12288 total bits
    for (int i = 0; i < 1024; i++) {
        bufs[i][0]= 0; //placeholder assignment
    }
    **/

    //We already have lines 1-3 from previous KeyGen part.

    //******************************STEP 4******************************
    //s′ = βμ(Rl×1q ; r)

    Vec<NTL::ZZ_pX> s_prime;
    s_prime.SetLength(l);

    k = 0;
    for (i = 0; i < l; i++) {
        for (j = 0; j < n; j++) {
            //we need to implement s[i][j] = HammingWeight(bufk) − HammingWeight(bufk+1) mod q
            ZZ_p temp = ZZ_p(buf_splitted[k].count()-buf_splitted[k+1].count()); //check if count returns hamming weight of bitset
            SetCoeff(s_prime[i],j,temp); 
            k = k + 2;
        }
    }

    //******************************STEP 5******************************
    //b′ = ((As′ + h) mod q) >> (Eq − Ep) ∈ Rl×1p

    Vec<ZZ_pX> b_vec_prime;
    b_vec_prime.SetLength(l);

    //this one is not transposed
    for (i = 0; i < l; i++) {
        for (j = 0; j < l; j++) {
            b_vec_prime[i] = matrix_A_NTL[i][j]*s_prime[i] + h_vec[i];
        }
    }

    for (i = 0; i < l; i++) {
        for (j = 0; j < l; j++) {
            //b_vec[i] >> (eq - EP);
            unsigned int temp3;
            conv(temp3, b_vec_prime[i][j]);
            temp3 >>= (eq - EP);
            b_vec_prime[i][j] = temp3;
        }
    }


    //******************************STEP 6******************************
    //v′ = b^T (s′ mod p) ∈ Rp

    ZZ_p::init(ZZ(p));
    ZZ_pX v_prime;

    for (i = 0; i < l; i++) {
        v_prime = b_vec_prime[i]*s_prime[i];
    }

    //******************************STEP 7******************************
    //c_m = (v′ + h1 − 2^(Ep−1)m mod p) >> (Ep − ET ) ∈ RT
    //c_m polynomial coefficients must be less than T = 8 for lightSABER, table 2 page 12

    ZZ_p::init(ZZ(2));
    ZZ_pX m; //this is the message m
    random(m, n);

    ZZ_p::init(ZZ(p));
    ZZ_pX temp_poly;
    temp_poly.SetLength(n);

    for (i = 0; i < n; i++) {
        SetCoeff(temp_poly, i, ZZ_p(pow(2,EP-1)*m[i]));
    }

    ZZ_pX c_m;
    //c_m = (v_prime + h1 - temp_poly) >> (EP - ET);
    c_m = (v_prime + h1 - temp_poly);

    for(size_t i = 0 ; i < n ; i++){
        unsigned int temp3;
        conv(temp3, c_m[i]);
        temp3 >>= (EP-ET);
        c_m[i] = temp3;
        //temp2[i] /= 512; // >> EP-1
    }

    //this is only adjusting the first 256 values to below T, all values must be below T
    for(size_t i = 0 ; i < n ; i++){
        long temp3;
        conv(temp3, c_m[i]);
        temp3 = temp3 % T;
        c_m[i] = temp3;
        //temp2[i] /= 512; // >> EP-1
    }

    //******************************STEP 8******************************
    //return c := (cm, b′)

    //****************************DECRYPTION****************************
    //s - Vector of L polynomials From Rq (same as Rp but with q as the p)
    //Rp set of some polynomial with coefficents less than P. P is a the modulus (is infinitly large)
    // b' Vector of polynomials in RP
    //s % p - 
    //R2 set of polynomials coefficents are less than 2 for example, 1 + x^2 + x^100     is in R_2

    /* Line 1 */
    
    // b' s and c come from other steps 
    //TODO - something in Decryption is modifying matrix_A output
    
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
    power(pow3,2,eq - EP - 1);
    
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

    Vec<ZZ_pX> s1; //had to change this to s1 because s is still a matrix
    Vec<ZZ_pX> b;
    
    s1.SetLength(l);
    b.SetLength(l);

    for(size_t i = 0 ; i < s1.length() ; i++){
        random(s1[i], degree);
        random(b[i], degree);
    }

    //b'T (s % p)
    for(size_t i = 0; i < l ; i++){
        v+= s1[i] * b[i];
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

    //std::cout << temp2 <<'\n';
    for(size_t i = 0 ; i < degree ; i++){
        unsigned int temp3;
        conv(temp3, temp2[i]);
        temp3 >>= (EP-1);
        temp2[i] = temp3;
        //temp2[i] /= 512; // >> EP-1
    }
    //std::cout << temp2 <<'\n';
    

    //******************************TESTS******************************
    //test seed generation from line 1 or 3
    cout << "Testing Seed Generation:" << endl;
    cout << "[";
    for(int i = 0; i < seed_A.size(); i++) {
        cout << seed_A[i];
    }
    cout << "]" << endl << endl;
    
    //test shake128 from line 2
    cout << "Testing SHAKE128 Output:" << endl;
    std::cout << "Buf: ";
    for (unsigned char c : buf){
            printf("%x",c);        
            }
    std::cout << std::endl;
    cout << endl;

    //test matrix generation from line 2 (per matrix element)
    cout << "Testing Matrix A Generation:" << endl;
    for(int i_1 = 0; i_1 < l; i_1++) {
        for(int i_2 = 0; i_2 < l; i_2++) {
            cout << matrix_A[i_1][i_2] << " ";
        } 
        cout << endl;
    }
    cout << endl;

    //test matrix generation from line 2 (per polynomial)
    for(int i_1 = 0; i_1 < l; i_1++) {
        for(int i_2 = 0; i_2 < l; i_2++) {
            for(int j = 0; j < n; j++) {
                cout << matrix_A[i_1][i_2][j] << " ";
            }
            cout << endl;
        } 
        cout << endl;
    }

    /**
    //testing step 4
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            for (k = 0; k < l; k++) {
                cout << s[i][j][k] << ' ';
            }
        }
        cout << endl;
    }
    **/

   /**
   //test results for step 5
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
    **/
    
    return 0;
}

//terminal arguments for compiling on Linux:
//g++ -I/usr/local/include SABER.cpp -o SABER.out -L/usr/local/lib -lcryptopp -lntl -lpthread
//./SABER.out 