#include <iostream>
#include <math.h>
#include <cmath>
#define PI 3.1415926
using namespace std;
int main()
{
  cout << "begining" << endl;
  // 6个晶粒，网格为100 * 100 
  const int N = 6;
  const int Nx = 100;
  const int Ny = 100;
  const int deltaX = 1;
  double phi[N][Nx][Ny]; // 当前时刻的phi
  double phi_b[N][Nx][Ny]; // 下一时刻的phi_b 

  // 初始化每个晶粒在不同网格处的phi值 
  int n, i, j;
  for (n = 0; n < N; n++)
  {
    for (i = 0; i < Nx; i++)
    {
      for (j = 0; j < Ny; j++)
      {
        phi[n][i][j] = phi_b[n][i][j] = 0; // 将晶粒在所有网格的phi都设置为0；

        // 将6个晶粒进行初始化，即特定位置的phi为1，其他位置都是0
        // 中间的三个进行初始化 
        if (n == 0 || n == 1 || n == 2) {
          if (j <= Ny / 2 + 3 && j >= Ny / 2 - 3){
            if (n == 0 && (i <= Nx / 6 + 3 && i >= Nx / 6 - 3)) {
              phi[0][i][j] = phi_b[0][i][j] = 1; // 第一个 
            }
            else if (n == 1 && (i <= (Nx / 6) * 3 + 3 && i >= (Nx / 6) * 3 - 3)) {
              phi[1][i][j] = phi_b[1][i][j] = 1; // 第二个 
            }
            else if (n == 2 && (i <= (Nx / 6) * 5 + 3 && i >= (Nx / 6) * 5 - 3)) {
              phi[2][i][j] = phi[2][i][j] = 1; // 第三个 
            }
          }
        }
        else if (n == 3 || n == 4) {
          // 上下的两个晶粒初始化
          if (j <= 3 || j >= Ny - 3) {
            if (n == 3 && (i <= Nx / 3 + 3 && i >= Nx / 3 - 3)) {
              phi[3][i][j] = phi_b[3][i][j] = 1; // 第四个 
            }
            else if (n == 4 && (i <= (Nx / 3) * 2 + 3 && i >= (Nx / 3) * 2 - 3)) {
              phi[4][i][j] = phi_b[4][i][j] = 1; // 第五个 
            }
          }
        }
        // 四个角的晶粒
        else if (n == 5){
          if ((i <= 3 && j <= 3) || (i >= Nx - 3 && j <= 3) || (i >= Nx - 3 && j >= Ny - 3) || (i <= 3 && j >= Ny - 3)) {
            phi[5][i][j] = phi_b[5][i][j] = 1; // 第六个 
          }
        }
      }
    }
  }



  // 测试初始化使用  
  // 通过以下程序可以测试出在初始化状态下每个晶粒的phi不为0的位置。 
  //for (int i = 0; i < nx; i++)
  //{
  //  for (int j = 0; j < ny; j++)
  //  {
  //    if (phi[2][i][j] == 1) {
  //      cout << i << '\t' << j << endl;
  //    }
  //  }
  //}

  //  cout << phi[3][62][198] << phi[3][70][1] << phi[3][67][196] << phi[3][68][2] << endl; // 测试初始化使用 
  //  cout << "上一个时刻第六个晶粒在原点处的phi值" << phi_b[5][0][0] << endl; // 测试phi_b的初始化使用 


    //double timeInterval = 0.1; // 记录时间步长
  double deltaT = 1; // 时间间
  int allTime = 100; // 晶粒生长的总时间
  double garma = 0.0208; // 单位 J/m2
  double thigma = 7 * 0.5 *  pow(10, -6); // thigma为7 * delta x， 而delta x为0.5um，这里我们换算为M做单位 
  double W = 4 * garma / thigma;
  double a = (2 / PI) * pow(2 * thigma * garma, 0.5);
  double deltaE = 0.09; // 单位Mpa
  double Qb = 110; // 单位 j/mol 
  double R = 8.314 * pow(10, -3); // 单位 kj/(K*mol) 
  double T = 800; // 这里的T为开尔文温度 
  double M = (0.139 / 800) * exp(-Qb / (R * T)); // M 0.139的单位是m2，是M0的表示 
  //cout << garma << endl << thigma << endl << W << endl << a << endl << deltaE << endl << Qb << endl << R << endl << T << endl << M << endl;    

  
  //cout << n << endl << i << endl << j << endl;

  int time = 0;
  while (time < allTime) {
    int k = 0; // 用于标记和n临界的晶粒
    for (n = 0; n < 6; n++) {
      for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
          double temp = 0;
          double dif = 0;
          for (k = 0; k < 6; k++) {
            if (k == n) {
              continue;
            }
            dif = (( phi[k][i + 1][j] + phi[k][i - 1][j] + phi[k][i][j + 1] + phi[k][i][j - 1] - 4 * phi[k][i][j]) - (phi[n][i + 1][j] + phi[n][i - 1][j] + phi[n][i][j + 1] + phi[n][i][j - 1] - 4 * phi[n][i][j])) / pow(deltaX, 2);
            temp += (M / N)*(W * (phi[k][i][j] - phi[n][i][j]) + 0.5 * pow(a, 2) * (dif) - (8 / PI) * pow(phi[n][i][j] * phi[k][i][j], 0.5) * deltaE);
          }
          phi_b[n][i][j] = temp * deltaT + phi[n][i][j];
        }
      }
    }
    time += deltaT;
  }


  cout << "ending" << endl;
  return 0;
}
