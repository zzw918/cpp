#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include <cstring>
#include <direct.h>
#include <windows.h>
#include "phi.h"
#define PI 3.1415926
using namespace std; 

extern const int N = 6; // add extern so phi.h can have access
extern const int Nx = 100; // Nx is pow(3, 0.5) times of Ny
extern const int Ny = 100; 

int ti(int); // deal with x-axis periodic boundry
int tj(int); // deal with y-axis periodic boundry

// difine the output function
int output(double phi_b[][Nx][Ny], int N, int Nx, int Ny, int fileNum);
int init(double phi[][Nx][Ny], double[][Nx][Ny], int N, int Nx, int Ny);
int foo(double phi[][Nx][Ny], double phi_b[][Nx][Ny], int n, int i, int j);

int main()
{
    cout << "begin to calculate." << endl;

    // define the phi array
    double phi[N][Nx][Ny];   // current phi
    double phi_b[N][Nx][Ny]; // next time phi

    // init the phi value
    init(phi, phi_b, N, Nx, Ny);
    output(phi_b, N, Nx, Ny, 100);

    // set the interval time and the whole time
    double deltaT = 0.01;                     // timeInterval
    double allTime = 5.0;                  // the whole time to grow

    double garma = 0.208;                 //  J/m2
    double deltaX = 0.5e-6;
    double Qb = 110e3;                               //  j/mol
    double R = 8.314;                              //  j/(K*mol)
    double T = 800.0;                                //  T is Kelvin's temperature
    double thigma = 7 * deltaX;           // delta x is 0.5um，use the 'm'
    double W = 4 * garma / thigma;
    double a = (2 / PI) * pow(2 * thigma * garma, 0.5);
    double M = (0.139 / T) * exp(-Qb / (R * T)) * PI * PI / (8 * thigma);
    cout << M << " " << W << " " << a << " ";

    double curTime = 0.0;
    int aid = 0;
    int number = 0;

    int files = 6;

    while (curTime < allTime)
    {   
        // every specified times to description. 
        if ((aid % (int(allTime / deltaT) / files)) == 0) {
            cout << double(curTime / allTime) * 100 << "% has calculated..." << endl;
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
                    for (k = 0; k < 6; k++)
                    {
                        dif = ((phi[k][ti(i + 1)][j] + phi[k][ti(i - 1)][j] + phi[k][i][tj(j + 1)] + phi[k][i][tj(j - 1)] - 4 * phi[k][i][j]) 
                        - (phi[n][ti(i + 1)][j] + phi[n][ti(i - 1)][j] + phi[n][i][tj(j + 1)] + phi[n][i][tj(j - 1)] - 4 * phi[n][i][j])) / pow(deltaX, 2);
                        temp += (2 * M / N) * (W * (phi[k][i][j] - phi[n][i][j]) + 0.5 * pow(a, 2) * dif);
                    }
                    phi_b[n][i][j] = -temp * deltaT + phi[n][i][j];
                    if (phi_b[n][i][j] < 10e-5) {
                        phi_b[n][i][j] = 0;
                    }
                }
            }
        }

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

        curTime += deltaT;
        if (aid % ((int)(allTime / deltaT) / files) == 0)
        {
            // 输出数据文件
            output(phi_b, N, Nx, Ny, number);
            number++;
        }
    }

    // output the file
    // output(phi_b, N, Nx, Ny, 0);

    cout << "calculate end!" << endl << endl;
    system("pause");
    return 0;
}



// 定义初始化函数
int init(double phi[][Nx][Ny], double phi_b[][Nx][Ny], int N, int Nx, int Ny) {
    int round = 2;
    for (int n = 0; n < N; n++)
    {
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                phi[n][i][j] = phi_b[n][i][j] = 0.165; // 0

                // 将6个晶粒进行初始化，即特定位置的phi为1，其他位置都是0
                // 中间的三个进行初始化
                if (n == 0 || n == 1 || n == 2)
                {
                    if (j <= Ny / 2.0 + round && j >= Ny / 2.0 - round)
                    {
                        // 注意：这里取6.0是因为int型相除得到的是int型，而只要分子或者分母有一个是double型，结果就是double型
                        if (n == 0 && (i <= Nx / 6.0 + round && i >= Nx / 6.0 - round))
                        {
                            phi[0][i][j] = phi_b[0][i][j] = 1.0; // 第一个
                            foo(phi, phi_b, n, i, j);
                        }
                        else if (n == 1 && (i <= (Nx / 6.0) * 3 + round && i >= (Nx / 6.0) * 3 - round))
                        {   
                            phi[1][i][j] = phi_b[1][i][j] = 1.0; // 第二个
                            foo(phi, phi_b, n, i, j);
                        }
                        else if (n == 2 && (i <= (Nx / 6.0) * 5 + round && i >= (Nx / 6.0) * 5 - round))
                        {
                            phi[2][i][j] = phi_b[2][i][j] = 1.0; // 第三个
                            foo(phi, phi_b, n, i, j);
                        }
                    }
                }
                else if (n == 3 || n == 4)
                {
                    // top and floor
                    if (j <= round || j >= Ny - round)
                    {
                        if (n == 3 && (i <= Nx / 3.0 + round && i >= Nx / 3.0 - round))
                        {
                            phi[3][i][j] = phi_b[3][i][j] = 1.0; // 第四个
                            foo(phi, phi_b, n, i, j);
                        }
                        else if (n == 4 && (i <= (Nx / 3.0) * 2 + round && i >= (Nx / 3.0) * 2 - round))
                        {
                            phi[4][i][j] = phi_b[4][i][j] = 1.0; //
                            foo(phi, phi_b, n, i, j);
                        }
                    }
                }

                // 四个角的晶粒
                else if (n == 5)
                {
                    if ((i <= round && j <= round) || (i >= Nx - round && j <= round) || (i >= Nx - round && j >= Ny - round) || (i <= round && j >= Ny - round))
                    {
                        phi[5][i][j] = phi_b[5][i][j] = 1.0; // 第六个
                        foo(phi, phi_b, n, i, j);
                    }
                }
            }
        }
    }

    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            double sum = 0.0;
            for (int n = 0; n < N; n++)
            {
                sum += phi[n][i][j];
            }
            if (sum > 0.0)
            {
                for (int n = 0; n < N; n++)
                {
                    double temp = phi[n][i][j];
                    phi[n][i][j] = temp / sum;
                }
            }
        }
    }
}

// 初始化相关函数
int foo(double phi[][Nx][Ny], double phi_b[][Nx][Ny], int n, int i, int j){
    for (int x = 0; x < N; x++) {
        if (x != n) {
            phi[x][i][j] = phi_b[x][i][j] = 0;
        }
    }
    return 0;
}


// 定义输出函数，方便重复调用，fileNum是用来生成不同的vtk文件的
int output(double phi_b[][Nx][Ny], int N, int Nx, int Ny, int fileNum) {

    double phi_o[Nx][Ny];
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            // 这里一定要是double，而不能是int，否则得到的结果则全部为整数，而显示出现问题！
            double sum = 0.0;
            for (int n = 0; n < N; n++)
            {
                sum = sum + pow(phi_b[n][i][j], 2);
            }
            phi_o[i][j] = sum;
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
    outfile << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << endl;
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
    return 0;
}