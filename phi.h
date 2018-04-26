#ifndef __PHI_H__
#define __PHI_H__
#include <iostream>
using namespace std;
extern const int N;
extern const int Nx;
extern const int Ny;

// deal with cyclical boundries
int ti(int i)
{
    // i will be out of bound
    if (i == -1) {
        return Nx - 2;
    } else if (i == Nx) {
        return 1;
    } else {
        return i;
    }
}

int tj(int j)
{
    // i will be out of bound
    if (j == -1) {
        return Ny - 2;
    } else if (j == Ny) {
        return 1;
    } else {
        return j;
    }
}

#endif