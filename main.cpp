#include <iostream>
#include <math.h>
#include <cmath>
#define PI 3.1415926
using namespace std;
int main()
{ 
    cout << "begining" << endl;
	// 6¸ö¾§Á££¬Íø¸ñÎª200 * 200 
	const int N = 6; 
	const int Nx = 100;
	const int Ny = 100;
	double phi[N][Nx][Ny]; // µ±Ç°Ê±¿ÌµÄphi
	double phi_b[N][Nx][Ny]; // ÉÏÒ»Ê±¿ÌµÄphi_b 
	
	// ³õÊ¼»¯Ã¿¸ö¾§Á£ÔÚ²»Í¬Íø¸ñ´¦µÄphiÖµ 
	int n, i, j;
	for (n = 0; n < N; n++)
	{
		for (i = 0; i < Nx; i++)
		{
			for (j = 0; j < Ny; j++)
			{ 
				phi[n][i][j] = phi_b[n][i][j] = 0; // ½«¾§Á£ÔÚËùÓÐÍø¸ñµÄphi¶¼ÉèÖÃÎª0£»
				
				// ½«6¸ö¾§Á£½øÐÐ³õÊ¼»¯£¬¼´ÌØ¶¨Î»ÖÃµÄphiÎª1£¬ÆäËûÎ»ÖÃ¶¼ÊÇ0
				// ÖÐ¼äµÄÈý¸ö½øÐÐ³õÊ¼»¯ 
				if (n == 0 || n == 1 || n == 2) {
					if (j <= Ny/2 + 3 && j >= Ny/2 - 3){
					    if (n == 0 && (i <= Nx/6 + 3 && i >= Nx/6 - 3)) {
					    	phi[0][i][j] = phi_b[0][i][j] = 1; // µÚÒ»¸ö 
					    }
					    else if (n == 1 && (i <= (Nx/6)*3 + 3 && i >= (Nx/6)*3 - 3 )) {
					    	phi[1][i][j] = phi_b[1][i][j] = 1; // µÚ¶þ¸ö 
					    }
					    else if (n == 2 && (i <= (Nx/6)*5 + 3 && i >= (Nx/6)*5 - 3 )) {
					    	phi[2][i][j] = phi[2][i][j] = 1; // µÚÈý¸ö 
				    	}
					} 
			    } 
				else if (n == 3 || n == 4) {
					// ÉÏÏÂµÄÁ½¸ö¾§Á£³õÊ¼»¯
				    if (j <= 3 || j >= Ny - 3) {
						if (n == 3 && (i <= Nx/3 + 3 && i >= Nx/3 - 3)) {
							phi[3][i][j] = phi_b[3][i][j] = 1; // µÚËÄ¸ö 
						}
						else if (n == 4 && (i <= (Nx/3)*2 + 3 && i >= (Nx/3)*2 - 3)) {
							phi[4][i][j] = phi_b[4][i][j] = 1; // µÚÎå¸ö 
						}
					} 
				}
				// ËÄ¸ö½ÇµÄ¾§Á£
				else if (n == 5){
					if ((i <= 3 && j <= 3) || (i >= Nx - 3 && j <= 3) || (i >= Nx - 3 && j >= Ny - 3) || (i <= 3 && j >= Ny - 3)) {
						phi[5][i][j] = phi_b[5][i][j] = 1; // µÚÁù¸ö 
					}
				}
			}
		}
	}
	
	

// ²âÊÔ³õÊ¼»¯Ê¹ÓÃ  
// Í¨¹ýÒÔÏÂ³ÌÐò¿ÉÒÔ²âÊÔ³öÔÚ³õÊ¼»¯×´Ì¬ÏÂÃ¿¸ö¾§Á£µÄphi²»Îª0µÄÎ»ÖÃ¡£ 
		for (int i = 0; i < Nx; i++)
		{
			for (int j = 0; j < Ny; j++)
			{
				if (phi[2][i][j] == 1) {
					cout << i << '\t' << j << endl;
				}
			}
		}

//  cout << phi[3][62][198] << phi[3][70][1] << phi[3][67][196] << phi[3][68][2] << endl; // ²âÊÔ³õÊ¼»¯Ê¹ÓÃ 
//	cout << "ÉÏÒ»¸öÊ±¿ÌµÚÁù¸ö¾§Á£ÔÚÔ­µã´¦µÄphiÖµ" << phi_b[5][0][0] << endl; // ²âÊÔphi_bµÄ³õÊ¼»¯Ê¹ÓÃ 

	
	int timeInterval = 0.1; // ¼ÇÂ¼Ê±¼ä²½³¤
	int allTime = 100; // ¾§Á£Éú³¤µÄ×ÜÊ±¼ä
	double garma = 0.0208; // µ¥Î» J/m2
	double thigma = 7 * 0.5 *  pow(10, -6); // thigmaÎª7 * delta x£¬ ¶ødelta xÎª0.5um£¬ÕâÀïÎÒÃÇ»»ËãÎªM×öµ¥Î» 
	double W = 4 * garma / thigma; 
	double a = (2 / PI) * pow(2 * thigma * garma, 0.5);
	double deltaE = 0.09; // µ¥Î»Mpa
	double Qb = 110; // µ¥Î» j/mol 
	double R = 8.314 * pow(10, -3); // µ¥Î» kj/(K*mol) 
	double T = 800; // ÕâÀïµÄTÎª¿ª¶ûÎÄÎÂ¶È 
	double M = (0.139 / 800) * exp(-Qb/(R * T)); // M 0.139µÄµ¥Î»ÊÇm2£¬ÊÇM0µÄ±íÊ¾ 
//	cout << garma << endl << thigma << endl << W << endl << a << endl << deltaE << endl << Qb << endl << R << endl << T << endl << M << endl;		 
	
	int time = 0; // 时间初始化

	

		
	cout << "ending" << endl;
	return 0;
}

