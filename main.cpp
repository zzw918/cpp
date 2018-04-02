#include <iostream>
#include <math.h>
#include <cmath>
#define PI 3.1415926
using namespace std;
int main()
{ 
    cout << "begining" << endl;
	// 6������������Ϊ200 * 200 
	const int N = 6; 
	const int Nx = 100;
	const int Ny = 100;
	double phi[N][Nx][Ny]; // ��ǰʱ�̵�phi
	double phi_b[N][Nx][Ny]; // ��һʱ�̵�phi_b 
	
	// ��ʼ��ÿ�������ڲ�ͬ���񴦵�phiֵ 
	int n, i, j;
	for (n = 0; n < N; n++)
	{
		for (i = 0; i < Nx; i++)
		{
			for (j = 0; j < Ny; j++)
			{ 
				phi[n][i][j] = phi_b[n][i][j] = 0; // �����������������phi������Ϊ0��
				
				// ��6���������г�ʼ�������ض�λ�õ�phiΪ1������λ�ö���0
				// �м���������г�ʼ�� 
				if (n == 0 || n == 1 || n == 2) {
					if (j <= Ny/2 + 3 && j >= Ny/2 - 3){
					    if (n == 0 && (i <= Nx/6 + 3 && i >= Nx/6 - 3)) {
					    	phi[0][i][j] = phi_b[0][i][j] = 1; // ��һ�� 
					    }
					    else if (n == 1 && (i <= (Nx/6)*3 + 3 && i >= (Nx/6)*3 - 3 )) {
					    	phi[1][i][j] = phi_b[1][i][j] = 1; // �ڶ��� 
					    }
					    else if (n == 2 && (i <= (Nx/6)*5 + 3 && i >= (Nx/6)*5 - 3 )) {
					    	phi[2][i][j] = phi[2][i][j] = 1; // ������ 
				    	}
					} 
			    } 
				else if (n == 3 || n == 4) {
					// ���µ�����������ʼ��
				    if (j <= 3 || j >= Ny - 3) {
						if (n == 3 && (i <= Nx/3 + 3 && i >= Nx/3 - 3)) {
							phi[3][i][j] = phi_b[3][i][j] = 1; // ���ĸ� 
						}
						else if (n == 4 && (i <= (Nx/3)*2 + 3 && i >= (Nx/3)*2 - 3)) {
							phi[4][i][j] = phi_b[4][i][j] = 1; // ����� 
						}
					} 
				}
				// �ĸ��ǵľ���
				else if (n == 5){
					if ((i <= 3 && j <= 3) || (i >= Nx - 3 && j <= 3) || (i >= Nx - 3 && j >= Ny - 3) || (i <= 3 && j >= Ny - 3)) {
						phi[5][i][j] = phi_b[5][i][j] = 1; // ������ 
					}
				}
			}
		}
	}
	
	

// ���Գ�ʼ��ʹ��  
// ͨ�����³�����Բ��Գ��ڳ�ʼ��״̬��ÿ��������phi��Ϊ0��λ�á� 
		for (int i = 0; i < Nx; i++)
		{
			for (int j = 0; j < Ny; j++)
			{
				if (phi[2][i][j] == 1) {
					cout << i << '\t' << j << endl;
				}
			}
		}

//  cout << phi[3][62][198] << phi[3][70][1] << phi[3][67][196] << phi[3][68][2] << endl; // ���Գ�ʼ��ʹ�� 
//	cout << "��һ��ʱ�̵�����������ԭ�㴦��phiֵ" << phi_b[5][0][0] << endl; // ����phi_b�ĳ�ʼ��ʹ�� 

	
	int timeInterval = 0.1; // ��¼ʱ�䲽��
	int allTime = 100; // ������������ʱ��
	double garma = 0.0208; // ��λ J/m2
	double thigma = 7 * 0.5 *  pow(10, -6); // thigmaΪ7 * delta x�� ��delta xΪ0.5um���������ǻ���ΪM����λ 
	double W = 4 * garma / thigma; 
	double a = (2 / PI) * pow(2 * thigma * garma, 0.5);
	double deltaE = 0.09; // ��λMpa
	double Qb = 110; // ��λ j/mol 
	double R = 8.314 * pow(10, -3); // ��λ kj/(K*mol) 
	double T = 800; // �����TΪ�������¶� 
	double M = (0.139 / 800) * exp(-Qb/(R * T)); // M 0.139�ĵ�λ��m2����M0�ı�ʾ 
//	cout << garma << endl << thigma << endl << W << endl << a << endl << deltaE << endl << Qb << endl << R << endl << T << endl << M << endl;		 
	
	

		
	cout << "ending" << endl;
	return 0;
}

