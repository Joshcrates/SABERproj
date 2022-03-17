#include <cstdio>

extern int ENC(int argc, char** argv);
extern int KEY_GEN(int argc, char** argv);

int main(int argc, char** argv){
    printf("Saber test\n");
    
    KEY_GEN(argc,argv);
    ENC(argc,argv);

    return 0;

}