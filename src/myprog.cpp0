#include <stdio.h>
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
using CryptoPP::SHAKE128;


using namespace std;
 
 int main() {
     
     using namespace CryptoPP;

     printf("Hello\n\n");
     
     std::string msg = "Yoda said, Do or do not. There is not try.";
     std::string digest;

     SHAKE128 hash;
     hash.Update((const byte*)msg.data(), msg.size());
     digest.resize(hash.DigestSize());
     hash.Final((byte*)&digest[0]);

     std::cout << "Message: " << msg << std::endl;

     std::cout << "Digest: ";
     //StringSource(digest, true, new Redirector(encoder));
     std::cout << std::endl;

     return 0;
 }

//g++ -I/usr/local/include myprog.cpp -o myprog.out -L/usr/local/lib -lcryptopp
//./myprog.out 