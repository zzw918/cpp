#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
//i and j in this program are all representing the grid
using namespace std;
void VTK(double **PHI, int Nx, int Ny, int xt, double dx, double dy)
{
	//输出vtk图，*PHI：sum(etas^2*);Nx:x方向格子数；Ny：y方向格子数；dx：x方向格子间距；dy：y方向格子间隔
	char fname[256];
	FILE *fp = NULL;				   //*fp的字符为空 file在源文件中插入当前源文件名
	sprintf(fname, "grain%d.vtk", xt); //printf是字符串格式化命令，主要功能是把格式化的数据写入某个字符串中
	int nz = 1;
	int non;
	double x, y, z;
	non = Nx * Ny * nz;
	ofstream outfile;
	outfile.open(fname);
	outfile << "# vtk DataFile Version 2.0" << endl;
	outfile << fname << endl;
	outfile << "ASCII " << endl;
	outfile << "DATASET STRUCTURED_GRID" << endl;
	outfile << "DIMENSIONS " << Nx << " " << Ny << " " << nz << endl;
	outfile << "POINTS " << non << " float" << endl;
	for (int i = 1; i <= Nx; i++)
	{
		for (int j = 1; j <= Ny; j++)
		{
			x = (i - 1) * dx;
			y = (j - 1) * dy;
			z = 0.0;
			outfile << x << " " << y << " " << z << endl;
		}
	}
	outfile << "POINT_DATA " << non << endl;
	outfile << "SCALARS CON float 1" << endl;
	outfile << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < Nx; i++)
	{

		for (int k = 0; k < Ny; k++)
		{
			if (k == 0)
				outfile << endl; //输出文件到数据里
			outfile << PHI[i][k] << "\t";
		}
	}
	cout << "输出成功" << endl;
	outfile.close();
}

int main()
{
	int Nx, Ny, N, aa, bb, igrain, nxx, nyy, Nn; //ndx,ndy; //Nx，x方向的个点数；Ny，y方向的个点数，N，晶粒数
	Nx = Ny = 128;								 //x和y方向格点数
	N = 4;
	Nn = 10000; //晶粒总数
	aa = 1, bb = 1;
	double rho0, delta_E, delta_x, Pi, TMAX, TIME, delta_t, gama, delta, W, M, T, Qb, M0, a, R, Mphi;
	vector<double> pvec(4, rho0);
	vector<double> rhoi(4, 0);
	double rhoc, step, delta_epsilon, epsilon, sigma_s, Qa, k2, A1, A2, delta_n, U, c, A, B, k1, sigama_c, tau, ngb, delta_y, d;
	double arr_p[1000][2] = {0};
	//    double p[4];
	gama = 0.208;		 //界面能γ
	delta_x = 0.5e-6;	//格点间距
	delta = 7 * delta_x; //界面厚度 δ
	delta_t = 0.01;		 //时间间隔
	TMAX = 30;			 //最长计算时间
	TIME = TMAX / delta_t;
	Pi = 3.1415926; //π
	W = 4 * gama / delta;
	T = 800; //温度
	M0 = 0.139;
	Qb = 110e3; //晶界扩散激活能
	R = 8.314;
	M = M0 / T * exp(-Qb / (R * T));  //晶界迁移率
	Mphi = Pi * Pi * M / (8 * delta); //场迁移率
	cout << "M=" << M << " "
		 << "Mphi=" << Mphi << endl;
	a = 2 * (sqrt(2 * delta * gama)) / Pi; //梯度能系数
	cout << "a=" << a << endl;
	delta_E = 0.0; //晶粒间储存能差值，

	rho0 = 10e9;   //初始位错密度
	k1 = 3.71e8;
	A = 0.5;		  //计算k2时用的系数
	A1 = 2.0e44;	  //计算sigma_s的系数
	A2 = 7.6;		  //计算sigma_s的系数
	U = 4.21e10;	  //切变模量 ，计算k2时使用
	B = 2.5e-10;	  //伯氏矢量
	epsilon = 2e-3;   //应变率
	Qa = 2.75e5;	  //激活能，计算sigma_s的系数
	d = 1;			  //计算形核率 delta_n时用的指数
	c = 5.0e25;		  //计算形核率 delta_n时用的系数
	sigama_c = 4.0e7; //计算rhco
	delta_epsilon = epsilon * delta_t;
	cout << "delta_epsilon=" << delta_epsilon << endl; //变形速率

	sigma_s = pow(((double)A1 * exp((double)Qa / R / T)), (double)1 / A2); //稳态应力
	cout << "sigma_s=" << sigma_s << endl;
	k2 = A * U * B * k1 / sigma_s; //动态回复的系数，
	cout << "k2=" << k2 << endl;
	rhoc = pow((sigama_c / A / U / B), 2); //临界位错密度
	cout << "rho=" << rhoc << endl;
	delta_n = (c * pow(delta_epsilon, d) * exp(-Qa / R / T)); //形核率
	cout << "delta_n=" << delta_n << endl;
	step = (delta_n * delta_t * ngb * delta_x * delta_y) / delta; // 步数

	//构造etas[][][]三维矩阵储存ith的φi

	double ***etas = new double **[N];
	for (int i = 0; i < N; i++)
	{
		etas[i] = new double *[Nx];
		for (int j = 0; j < Nx; j++)
		{
			etas[i][j] = new double[Ny];
		}
	}
	double dx, dy; //画坐标
	dx = dy = 0.5;
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int igrain = 0; igrain < 4; igrain++)
			{
				etas[igrain][i][j] = 1.0 / 4.0;
			}
		}
	}
	double ***etas1 = new double **[N];
	for (int i = 0; i < N; i++)
	{
		etas1[i] = new double *[Nx];
		for (int j = 0; j < Nx; j++)
		{
			etas1[i][j] = new double[Ny];
		}
	}
	double ***delta_etas = new double **[N];
	for (int i = 0; i < N; i++)
	{
		delta_etas[i] = new double *[Nx];
		for (int j = 0; j < Nx; j++)
		{
			delta_etas[i][j] = new double[Ny];
		}
	}
	for (int i = 0; i < Nx / aa; i++)
	{
		for (int j = 0; j < Ny / bb; j++)
		{
			for (int igrain = 0; igrain < 4; igrain++)
			{
				if (igrain == 0 && pow((i - (Nx - 1) / aa / 2), 2) + pow((j - (Ny - 1) / bb / 2), 2) <= pow(5, 2))
				{
					etas[igrain][i][j] = 1;
					for (int m = 0; m < N; m++)
					{
						if (m != igrain)
						{
							etas[m][i][j] = 0;
						}
					}
				}
				if (igrain == 1 && pow(i / aa, 2) + pow(j - ((Ny - 1) / bb / 4), 2) <= pow(5, 2))
				{
					etas[igrain][i][j] = 1;
					for (int m = 0; m < N; m++)
					{
						if (m != igrain)
						{
							etas[m][i][j] = 0;
						}
					}
				}
				if (igrain == 2 && pow(i / aa, 2) + pow((j - 3 * (Ny - 1) / bb / 4), 2) <= pow(5, 2))
				{
					etas[igrain][i][j] = 1;
					for (int m = 0; m < N; m++)
					{
						if (m != igrain)
						{
							etas[m][i][j] = 0;
						}
					}
				}
				if (igrain == 3 && pow((i - (Nx - 1) / aa / 2), 2) + pow(j / bb, 2) <= pow(5, 2))
				{
					etas[igrain][i][j] = 1;
					for (int m = 0; m < N; m++)
					{
						if (m != igrain)
						{
							etas[m][i][j] = 0;
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int g = 0; g < N; g++)
			{
				int a, b;
				a = i;
				b = j;
				if (i >= (Nx / aa - 1))
				{
					a = i - (Nx / aa - 1);
				}
				if (j >= (Ny / aa - 1))
				{
					b = j - (Ny / aa - 1);
				}
				etas[g][i][j] = etas[g][a][b];
			}
		}
	}
	double **PHI = new double *[Nx];
	for (int i = 0; i < Nx; i++)
	{
		PHI[i] = new double[Ny];
	}
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			double sum = 0;
			for (int n = 0; n < N; n++)
			{
				sum = sum + etas[n][i][j] * etas[n][i][j];
			}
			PHI[i][j] = sum;
		}
	}
	VTK(PHI, Nx, Ny, 0, dx, dy);
	cout << "进行相场计算,计算时长TIME=" << TMAX << endl;
	for (int Time = 1; Time <= TIME; Time++)
	{
		for (int x = 0; x < Nx; x++)
		{
			for (int y = 0; y < Ny; y++)
			{
				for (int m = 0; m < N; m++) //边界条件
				{
					int a, b, c, d;
					a = x + 1;
					b = x - 1;
					c = y + 1;
					d = y - 1;
					if (a == Nx)
					{
						a = 1;
					}
					if (b < 0)
					{
						b = Nx - 2; //
					}
					if (c == Ny)
					{
						c = 1;
					}
					if (d < 0)
					{
						d = Ny - 2;
					}
					delta_etas[m][x][y] = (etas[m][a][y] + etas[m][b][y] + etas[m][x][c] + etas[m][x][d] - 4 * etas[m][x][y]) / (delta_x * delta_x);
				}
			}
		}
		for (int i = 0; i < N; i++)
		{
			for (int x = 0; x < Nx; x++)
			{
				for (int y = 0; y < Ny; y++)
				{
					double sum = 0.0;
					int n = 0;
					for (int j = 0; j < N; j++)
					{
						if (etas[j][x][y] != 0)
						{
							sum = sum + 2 * Mphi * (W * (etas[j][x][y] - etas[i][x][y]) + a * a * (delta_etas[j][x][y] - delta_etas[i][x][y]) / 2 - 8 * delta_E * sqrt(etas[i][x][y] * etas[j][x][y]) / Pi);
							n = n + 1;
						}
					}
					etas1[i][x][y] = -sum * delta_t / n + etas[i][x][y];
				}
			}
		}
	}

	//再结晶
	cout << "进行再结晶计算,计算时长TIME=" << TMAX << endl; //将下一时刻的rho值赋值给每一个晶粒
	for (int Time = 1; Time <= TIME; Time++)
	{
		for (int m = 0; m < pvec.size(); m++)
		{
			rhoi[m] = (k1 * sqrt(rho0) - k2 * rho0) * delta_epsilon + pvec[m];
			//	  cout<<p[0]<<p[1]<<p[2]<<p[3]<<endl;
			rhoi.push_back(rhoi[m]);
		}
		int ngb = 0, tx, ty, vi, vj;
		double sum = 0.0;
		for (int ii = 0; ii < N; ii++)
		{
			for (vi = 0; vi < Nx; vi++)
			{
				for (vj = 0; vj < Ny; vj++)
				{
					etas1[ii][vi][vj] = pow(etas1[ii][vi][vj], 2);
					sum = sum + etas1[ii][vi][vj];
				}
			}
		}
		for (int m = 0; m < pvec.size(); m++)
			if (rhoi[m] > rhoc && sum < 0.6) //判断再结晶条件
			{
				arr_p[ngb][0] = vi; //将满足再结晶的vi坐标赋值给arr_p[ngb][0]
				arr_p[ngb][1] = vj; //将满足再结晶的vj坐标赋值给arr_p[ngb][1]
				ngb++;				//ngb为一共是几对横纵坐标
									//	        cout<<ngb<<" "<<vi<<" "<<vj<<endl;
			}
		int nn = 1, num;				 //nn是满足step大于整数时形核
		srand(time(0));					 //随机的从p数组中选取成对的x，y
		for (int t = 0; t < 100000; t++) //
		{
			if (step > nn)
			{
				num = rand() % (ngb);
				// arr_p[num][0];
				// arr_p[num][1];
				// 这里写的还是非常棒的，我之前都没有想到，当然，用vector也是可以实现的。
				int tx = arr_p[num][0], ty = arr_p[num][1];
				for (int iix = 0; iix < Nx / aa; iix++)
				{
					for (int iiy = 0; iiy < Ny / bb; iiy++)
					{
						for (int igrain = 4; igrain < Nn; igrain++)
							if (pow((iix - tx), 2) + pow((iiy - ty), 2) <= pow(5 * delta_x, 2))
							{
								etas[igrain][iix][iiy] = 1;
								//		    	        	for(int m_vec=4;m_vec<Nx;m_vec++)
								//		    	        	{

								//		    	        	}
								igrain++;
								pvec.push_back(rho0);
							}
					}
				}
				nn++;
			}
		}
		for (int i = 0; i < N; i++)
		{
			for (int x = 0; x < Nx; x++)
			{
				for (int y = 0; y < Ny; y++)
				{
					double sum = 0.0;
					int n = 0;
					for (int j = 0; j < N; j++)
					{
						if (etas[j][x][y] != 0)
						{
							sum = sum + 2 * Mphi * (W * (etas[j][x][y] - etas[i][x][y]) + a * a * (delta_etas[j][x][y] - delta_etas[i][x][y]) / 2 - 8 * delta_E * sqrt(etas[i][x][y] * etas[j][x][y]) / Pi);
							n = n + 1;
						}
					}
					etas1[i][x][y] = -sum * delta_t / n + etas[i][x][y];
				}
			}
		}
		for (int i = 0; i < Nn; i++)
		{
			for (int x = 0; x < Nx; x++)
			{
				for (int y = 0; y < Ny; y++)
				{
					if (etas1[i][x][y] <= 1.0e-5)
					{
						etas1[i][x][y] = 0.0;
					}
				}
			}
		}
		for (int x = 0; x < Nx; x++)
		{
			for (int y = 0; y < Ny; y++)
			{
				double sum = 0.0;
				for (int i = 0; i < Nn; i++)
				{
					sum = sum + etas1[i][x][y];
				}
				for (int i = 0; i < Nn; i++)
				{
					etas[i][x][y] = etas1[i][x][y] / sum;
				}
			}
		}
		double tttt, temp;
		tttt = Time / TIME * 100;
		temp = fmod(tttt, 10);
		if (temp == 0)
		{
			if (tttt == 100)
			{
				cout << ".....计算完成，输出VTK文件.....";
			}
			else
			{
				cout << "计算中，已完成"
					 << "\t" << tttt << ".0%" << endl;
			}
			char fname[256];
			int a;
			a = tttt / 10;
			double **PHI = new double *[Nx];
			for (int i = 0; i < Nx; i++)
			{
				PHI[i] = new double[Ny];
			}
			for (int i = 0; i < Nx; i++)
			{
				for (int j = 0; j < Ny; j++)
				{
					double sum = 0;
					for (int n = 0; n < Nn; n++)
					{
						sum = sum + etas[n][i][j] * etas[n][i][j];
					}
					PHI[i][j] = sum;
				}
			}
			VTK(PHI, Nx, Ny, a, dx, dy);
		}
	}
	return 0;
}