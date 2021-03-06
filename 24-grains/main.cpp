#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <direct.h>
#include <windows.h>
#define PI 3.1415926
using namespace std; 

const int N = 25; // the total number of grains
// that Nx is 120 and Ny is 70 is fit
const int Nx = 256; // the x-axis grid numbers
const int Ny = 150; // the y-axis grid numbers

int ti(int); // deal with x-axis periodic boundry
int tj(int); // deal with y-axis periodic boundry
void output(double phi_b[][Nx][Ny], int N, int Nx, int Ny, int fileNum); // difine the output function
void init(double phi[][Nx][Ny], double phi_b[][Nx][Ny], int N, int Nx, int Ny); // define the initiated function
void init_zero(double phi[][Nx][Ny], double phi_b[][Nx][Ny], int n, int i, int j); // init the area where is not belong to the core

int main()
{
    cout << "begin to calculate." << endl;

    // define the phi array, use the static array, avoid the overflow
    static double phi[N][Nx][Ny];   // current phi
    static double phi_b[N][Nx][Ny]; // next time phi
    // init the phi value
    init(phi, phi_b, N, Nx, Ny);

    // TEST, test the data, look if it is right
    // for (int n = 0; n < N; n++) {
    //     for (int i = 0; i < Nx; i++) {
    //         for (int j = 0; j < Ny; j++) {
    //             if (phi[n][i][j] == 1) {
    //                 cout << n << " " << i << " " << j << endl;
    //             }
    //         }
    //     }
    // }

    // output the initiated file here
    output(phi_b, N, Nx, Ny, 0);



    // set the interval time and the whole time
    // for 24 grains, that daletT is 0.01 and allTime is 4.0 is fit
    double deltaT = 0.01,                 // timeInterval
           allTime = 30.0;                  // the whole time to grow

    double garma = 0.208,                  //  J/m2
           deltaX = 0.5e-6,
           Qb = 110e3,                     //  j/mol
           R = 8.314,                      //  j/(K*mol)
           T = 800.0,                      //  T is Kelvin's temperature
           thigma = 7 * deltaX,            // delta x is 0.5um，use the 'm'
           W = 4 * garma / thigma,
           a = (2 / PI) * pow(2 * thigma * garma, 0.5),
           M = (0.139 / T) * exp(-Qb / (R * T)) * PI * PI / (8 * thigma);

    double curTime = 0.0;
    int aid = 0; // to help show output
    int number = 1; // filename number mark

    int files = 8; // file numbers

    while (curTime < allTime)
    {   
        // every specified times to description. 
        if ((aid % (int(allTime / deltaT) / files)) == 0) {
            cout << fixed << setprecision(0) << double(curTime / allTime) * 100 << "% has been calculated..." << endl;
        }

        aid++;

        int k = 0; // mark the grain which is closed to the n grain
        for (int n = 0; n < N; n++)
        {   
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    double temp = 0.0;
                    double dif = 0.0;
                    double E = 0.0;
                    for (k = 0; k < N; k++)
                    {
                        if (n == N - 1 && k != N - 1) {
                            E = -0.09e6;
                        }
                        if (n != N - 1 && k == N - 1) {
                            E = +0.09e6;
                        }
                        dif = ((phi[k][ti(i + 1)][j] + phi[k][ti(i - 1)][j] + phi[k][i][tj(j + 1)] + phi[k][i][tj(j - 1)] - 4 * phi[k][i][j]) 
                        - (phi[n][ti(i + 1)][j] + phi[n][ti(i - 1)][j] + phi[n][i][tj(j + 1)] + phi[n][i][tj(j - 1)] - 4 * phi[n][i][j])) / pow(deltaX, 2);
                        temp += (2 * M / N) * (W * (phi[k][i][j] - phi[n][i][j]) + 0.5 * pow(a, 2) * dif - 8/PI*pow(phi[n][i][j] * phi[k][i][j], 0.5) * E);
                    }
                    phi_b[n][i][j] = -temp * deltaT + phi[n][i][j];
                    if (phi_b[n][i][j] < 10e-5) {
                        phi_b[n][i][j] = 0;
                    }
                }
            }
        }

        // make sure that the sum of grains of one grid is 1
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                double sum = 0.0;
                for (int n = 0; n < N; n++)
                {
                    sum += phi_b[n][i][j];
                }
                if (sum > 0.0)
                {
                    for (int n = 0; n < N; n++)
                    {
                        double temp = phi_b[n][i][j];
                        phi_b[n][i][j] = temp / sum;
                    }
                }
            }
        }

        // assign the phi value from the phi_b value
        for (int n = 0; n < N; n++)
        {
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    phi[n][i][j] = phi_b[n][i][j];  
                }
            }
        }

        curTime += deltaT; // time added

        if (aid % ((int)(allTime / deltaT) / files) == 0)
        {
            // output every period data
            output(phi_b, N, Nx, Ny, number++);
        }
    }
    cout << "calculate end!" << endl << endl;
    system("pause");
    return 0;
}



// define the initiated function
void init(double phi[][Nx][Ny], double phi_b[][Nx][Ny], int N, int Nx, int Ny) {
    int round = 7;
    // first, make all the grain's phi be 1.0/24.0
    for (int n = 0; n < N; n++) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                phi[n][i][j] = phi_b[n][i][j] = 0;
                if (n == N - 1) {
                    phi[n][i][j] = phi_b[n][i][j] = 1;
                }
            }
        }
    }

    for (int n = 0; n < N; n++)
    {
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                // init the middle three grains
                switch (n) {
                    // new
                    // grains 1 - 6
                    case 0:
                        if ((i <= Nx / 12.0 + round && i >= Nx / 12.0 - round) && (j <= (Ny / 4.0) * 3 + round && j >= (Ny / 4.0) * 3 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the first one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;
                        
                    case 1:
                        if ((i <= (Nx / 12.0) * 3 + round && i >= (Nx / 12.0) * 3 - round) && (j <= (Ny / 4.0) * 3 + round && j >= (Ny / 4.0) * 3 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the second one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 2:
                        if ((i <= (Nx / 12.0) * 5 + round && i >= (Nx / 12.0) * 5 - round) && (j <= (Ny / 4.0) * 3 + round && j >= (Ny / 4.0) * 3 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the third one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 3:
                        if ((i <= (Nx / 12.0) * 7 + round && i >= (Nx / 12.0) * 7 - round) && (j <= (Ny / 4.0) * 3 + round && j >= (Ny / 4.0) * 3 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the third one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 4:
                        if ((i <= (Nx / 12.0) * 9 + round && i >= (Nx / 12.0) * 9 - round) && (j <= (Ny / 4.0) * 3 + round && j >= (Ny / 4.0) * 3 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the third one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 5:
                        if ((i <= (Nx / 12.0) * 11 + round && i >= (Nx / 12.0) * 11 - round) && (j <= (Ny / 4.0) * 3 + round && j >= (Ny / 4.0) * 3 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the third one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    // grains 7 - 11
                    case 6:
                        if ((i <= (Nx / 12.0) * 2 + round && i >= (Nx / 12.0) * 2 - round) && (j <= (Ny / 4.0) * 2 + round && j >= (Ny / 4.0) * 2 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the third one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 7:
                        if ((i <= (Nx / 12.0) * 4 + round && i >= (Nx / 12.0) * 4 - round) && (j <= (Ny / 4.0) * 2 + round && j >= (Ny / 4.0) * 2 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the third one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 8:
                        if ((i <= (Nx / 12.0) * 6 + round && i >= (Nx / 12.0) * 6 - round) && (j <= (Ny / 4.0) * 2 + round && j >= (Ny / 4.0) * 2 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the third one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 9:
                        if ((i <= (Nx / 12.0) * 8 + round && i >= (Nx / 12.0) * 8 - round) && (j <= (Ny / 4.0) * 2 + round && j >= (Ny / 4.0) * 2 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the third one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 10:
                        if ((i <= (Nx / 12.0) * 10 + round && i >= (Nx / 12.0) * 10 - round) && (j <= (Ny / 4.0) * 2 + round && j >= (Ny / 4.0) * 2 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the third one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    // grains 12 - 17
                    case 11:
                        if ((i <= Nx / 12.0 + round && i >= Nx / 12.0 - round) && (j <= (Ny / 4.0) * 1 + round && j >= (Ny / 4.0) * 1 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the first one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 12:
                        if ((i <= (Nx / 12.0) * 3 + round && i >= (Nx / 12.0) * 3 - round) && (j <= (Ny / 4.0) * 1 + round && j >= (Ny / 4.0) * 1 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the first one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 13:
                        if ((i <= (Nx / 12.0) * 5 + round && i >= (Nx / 12.0) * 5 - round) && (j <= (Ny / 4.0) * 1 + round && j >= (Ny / 4.0) * 1 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the first one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 14:
                        if ((i <= (Nx / 12.0) * 7 + round && i >= (Nx / 12.0) * 7 - round) && (j <= (Ny / 4.0) * 1 + round && j >= (Ny / 4.0) * 1 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the first one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 15:
                        if ((i <= (Nx / 12.0) * 9 + round && i >= (Nx / 12.0) * 9 - round) && (j <= (Ny / 4.0) * 1 + round && j >= (Ny / 4.0) * 1 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the first one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 16:
                        if ((i <= (Nx / 12.0) * 11 + round && i >= (Nx / 12.0) * 11 - round) && (j <= (Ny / 4.0) * 1 + round && j >= (Ny / 4.0) * 1 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the first one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    // grains 18 - 22
                    case 17:
                        if ((i <= Nx / 6.0 + round && i >= Nx / 6.0 - round) && (j <= round || j >= Ny - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the forth one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 18:
                        if ((i <= (Nx / 6.0) * 2 + round && i >= (Nx / 6.0) * 2 - round) && (j <= round || j >= Ny - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the forth one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 19:
                        if ((i <= (Nx / 6.0) * 3 + round && i >= (Nx / 6.0) * 3 - round) && (j <= round || j >= Ny - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the forth one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 20:
                        if ((i <= (Nx / 6.0) * 4 + round && i >= (Nx / 6.0) * 4 - round) && (j <= round || j >= Ny - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the forth one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    case 21:
                        if ((i <= (Nx / 6.0) * 5 + round && i >= (Nx / 6.0) * 5 - round) && (j <= round || j >= Ny - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the forth one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;
                    
                    // grain 23
                    case 22:
                        if ((i <= round || i >= Nx - round) && (j <= (Ny / 4.0) * 2 + round && j >= (Ny / 4.0) * 2 - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the third one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;

                    // grain 24
                    case 23:
                        if ((i <= round && j <= round) || (i >= Nx - round && j <= round) || (i >= Nx - round && j >= Ny - round) || (i <= round && j >= Ny - round))
                        {
                            phi[n][i][j] = phi_b[n][i][j] = 1.0; // the sixth one
                            init_zero(phi, phi_b, n, i, j);
                        }
                        break;
                }
            }
        }
    }
}

// init the area where is not belong to the core
void init_zero(double phi[][Nx][Ny], double phi_b[][Nx][Ny], int n, int i, int j){
    for (int x = 0; x < N; x++) {
        if (x != n) {
            phi[x][i][j] = phi_b[x][i][j] = 0.0;
        }
    }
}


// define the output function
void output(double phi_b[][Nx][Ny], int N, int Nx, int Ny, int fileNum) {

    double phi_o[Nx][Ny];
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            // HINT: here sum is "double" type not "int" type
            double sum = 0.0;
            for (int n = 0; n < N; n++)
            {
                sum = sum + pow(phi_b[n][i][j], 2);
            }
            phi_o[i][j] = sum;
            if (phi_o[i][j] > 1.0) {
                phi_o[i][j] = 1.0;
            }
        }
    }

    // TEST: here we use the data to test the result
    // for (int n = 0; n < N; n++)
    // {
    //     for (int i = 0; i < Nx; i++)
    //     {
    //         for (int j = 0; j < Ny; j++)
    //         {
    //             if (phi_b[n][i][j] == 1)
    //             {
    //                 cout << n << " " << i << " " << j << endl;
    //             }
    //         }
    //     }
    // }

    // output the vtk files
    // fstream outfile;
    fstream outfile;
    char filename[20] = "output";
    char num[20];
    char foldername[20] = "./outputfolder/";
    // _mkdir(foldername); // not completed!
    itoa(fileNum, num, 10);
    strcat(strcat(filename, num), ".vtk");
    outfile.open(strcat(foldername, filename), ios::out);
    if (!outfile)
    {
        cout << "open failed" << endl;
    }
    outfile << "# vtk DataFile Version 2.0" << endl;
    outfile << filename << endl;
    outfile << "ASCII " << endl;
    outfile << "DATASET STRUCTURED_GRID" << endl;
    outfile << "DIMENSIONS " << Ny << " " << Nx << " " << 1 << endl;
    outfile << "POINTS " << Nx * Ny * 1 << " float" << endl;
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            outfile << i << " " << j << " " << 1 << endl;
        }
    }

    outfile << "POINT_DATA " << Nx * Ny * 1 << endl;
    outfile << "SCALARS CON float 1" << endl;
    outfile << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            outfile << phi_o[i][j] << "\t";
        }
    }
    outfile.close();
}

int ti(int i)
{
    // i will be out of bound
    if (i == -1)
    {
        return Nx - 2;
    }
    else if (i == Nx)
    {
        return 1;
    }
    else
    {
        return i;
    }
}

int tj(int j)
{
    // j will be out of bound
    if (j == -1)
    {
        return Ny - 2;
    }
    else if (j == Ny)
    {
        return 1;
    }
    else
    {
        return j;
    }
}