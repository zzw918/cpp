#include <iostream>
#include <math.h>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <direct.h>
#include <windows.h>

#define PI 3.1415926

using namespace std;

const int N = 1;// the total number of grains
const int Nx = 50 * 2;// the x-axis grid numbers
const int Ny =50 * 2;// the y-axis grid numbers

int ti(int); // deal with x-axis periodic boundry
int tj(int); // deal with y-axis periodic boundry

void output(double phi_b[][Nx][Ny], int N, int Nx, int Ny, int fileNum); // difine the output function
void init(double phi[][Nx][Ny], double phi_b[][Nx][Ny], int N, int Nx, int Ny); // define the initiated function

int main()
{
cout << "begin to calculate." << endl;

    // define the phi array
    static double phi[N][Nx][Ny];   // current phi
    static double phi_b[N][Nx][Ny]; // next time phi

    // init the phi value
    init(phi, phi_b, N, Nx, Ny);

    // output the initiated file here
    output(phi_b, N, Nx, Ny, 0);

    // set the interval time and the whole time
    double deltaT = 0.01,                 // timeInterval
           allTime = 40.0;                  // the whole time to grow

    double garma = 0.208,                  //  J/m2
           deltaE = 0.09,
           deltaX = 0.5e-6,
           Qb = 110e3,                     //  j/mol
           R = 8.314,                      //  j/(K*mol)
           T = 800.0,                      //  T is Kelvin's temperature
           thigma = 7 * deltaX,            // delta x is 0.5um，use the'm'
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
                    phi_b[n][i][j] = ((-temp * deltaT)) + phi[n][i][j];

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
    double radius = 5.48;
    for (int n = 0; n < N; n++) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                phi[n][i][j] = phi_b[n][i][j] = 0;
                if (radius > sqrt( j*j + i*i )) {
                    phi_b[n][i][j] = phi[n][i][j] = 1;
                }   
            }
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
