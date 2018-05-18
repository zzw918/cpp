#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<fstream>
#include<cmath>

using namespace std;
void VTK(double **PHI, int Nx,int Ny,int i,double dx,double dy)
{
	//���vtkͼ��*PHI��sum(etas^2*);Nx:x�����������Ny��y�����������dx��x������Ӽ�ࣻdy��y������Ӽ��
	char fname[256];
	FILE *fp = NULL;   //*fp���ַ�Ϊ�� file��Դ�ļ��в��뵱ǰԴ�ļ���
	sprintf(fname, "grain%d.vtk",i);  //printf���ַ�����ʽ�������Ҫ�����ǰѸ�ʽ��������д��ĳ���ַ����� 
	int nz=1;
	int non;
	double x,y,z;
	non=Nx*Ny*nz;
	ofstream outfile;
	outfile.open(fname);
	outfile<<"# vtk DataFile Version 2.0" <<endl;
	outfile<<fname<<endl;
	outfile<<"ASCII "<<endl;
	outfile<<"DATASET STRUCTURED_GRID"<<endl;
	outfile<<"DIMENSIONS "<<Nx<<" "<<Ny<<" "<<nz<<endl;
	outfile<<"POINTS "<<non<<" float"<<endl;
	for(int i=1;i<=Nx;i++)
	{
		for(int j=1;j<=Ny;j++)
		{
			x=(i-1)*dx;
			y=(j-1)*dy;
			z=0.0;
			outfile<<x<<" "<<y<<" "<<z<<endl;
		}
	}
	outfile<<"POINT_DATA "<<non<<endl;
	outfile<<"SCALARS CON float 1"<<endl;
	outfile<<"LOOKUP_TABLE default"<<endl;
	for(int i=0;i<Nx;i++)
	{ 
		
		for(int k=0;k<Ny;k++)
		{
			if(k==0)
			outfile<<endl; //����ļ��������� 
			outfile<<PHI[i][k]<<"\t";
	    }
	}
	cout<<"����ɹ�"<<endl;;
	outfile.close();
}


int main()
{ 
	int Nx,Ny,N; //Nx��x����ĸ�������Ny��y����ĸ�������N�������� 
	Nx=Ny=128;//x��y�������� 
	N=7;//�������� 
	double delta_E,delta_x,Pi,TMAX,TIME,delta_t,gama,delta,W,M,T,Qb,M0,a,R,Mphi;
	gama=0.208;//�����ܦ� 
	delta_x=0.5e-6;//����� 
	delta=7*delta_x;//������ �� 
	delta_t=0.01;//ʱ���� 
	TMAX=30;//�����ʱ�� 
	TIME=TMAX/delta_t;
	Pi=3.1415926;//�� 
	W=4*gama/delta;
	T=800;//�¶� 
	M0=0.139;
	Qb=110e3;
	R=8.314;
	M=M0/T*exp(-Qb/(R*T));//����Ǩ���� 
	Mphi=Pi*Pi*M/(8*delta);//��Ǩ���� 
	cout<<M<<" "<<Mphi<<endl; 
	a=2*(sqrt(2*delta*gama))/Pi;//�ݶ���ϵ�� 
	delta_E=0.0;//�����䴢���ܲ�ֵ��
	//����etas[][][]��ά���󴢴�ith�Ħ�i 
	double ***etas=new double **[N];
	for(int i=0;i<N;i++)
	{
		etas[i]=new double *[Nx];
		for(int j=0;j<Nx;j++)
		{
			etas[i][j]=new double[Ny];
		}
	} 
	double dx,dy;
	dx=dy=0.5;
	double *gx=new double[Nx];
	double *gy=new double[Ny];
	for(int i=0;i<Nx;i++)
	{
		gx[i]=(i+1)*dx;
		gy[i]=(i+1)*dy;
	}
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Ny;j++)
		{
			for(int igrain=0;igrain<7;igrain++)
			{
				etas[igrain][i][j]=0.0;
			}
		}
	}
	double ***etas1=new double **[N];
	for(int i=0;i<N;i++)
	{
		etas1[i]=new double *[Nx];
		for(int j=0;j<Nx;j++)
		{
			etas1[i][j]=new double[Ny];
		}
	}
	double ***delta_etas=new double **[N];
	for(int i=0;i<N;i++)
	{   
		delta_etas[i]=new double *[Nx];
		for(int j=0;j<Nx;j++)
		{
			delta_etas[i][j]=new double[Ny];
		}
	}

	   for(int i=0;i<Nx;i++)
        {
	     for(int j=0;j<Ny;j++)
		 {
		   for(int igrain=0;igrain<7;igrain++)
			                             
		    {
		      if(igrain==0&&pow((i-(Nx-1)/2),2)+pow((j-(Ny-1)/2),2)<=pow(5,2)) 
			   {
			    etas[igrain][i][j]=1;
			   } 
			    if(igrain==1&&pow(i,2)+pow(j-((Ny-1)/4),2)<=pow(5,2)) 
			    {
                 etas[igrain][i][j]=1;	
			    }
                  if(igrain==2&&pow(i,2)+pow((j-3*(Ny-1)/4),2)<=pow(5,2)) 
                  {
	               etas[igrain][i][j]=1;
                  }
	              if(igrain==3&&pow((i-(Nx-1)/2),2)+pow(j,2)<=pow(5,2)) 
	              {
	              etas[igrain][i][j]=1;	
	              }
	              if(igrain==4&&pow((i-(Nx-1)/2),2)+pow((j-(Ny-1)),2)<=pow(5,2)) 
	              {
                   etas[igrain][i][j]=1;
	              }
                  if(igrain==5&& pow((i-(Nx-1)),2)+pow((j-(Ny-1)/4),2)<=pow(5,2)) 
                  {
	               etas[igrain][i][j]=1;
                  }
	              if(igrain==6&&pow((i-(Nx-1)),2)+pow(j-3*(Ny-1)/4,2)<=pow(5,2)) 
	              {
                   etas[igrain][i][j]=1;	
	              }  
		    }
		 }
        }
	cout<<"�����ೡ����,����ʱ��TIME="<<TMAX<<endl; 
	for(int time=1;time<=TIME;time++)
	{
		for(int x=0;x<Nx;x++)
		{
			for(int y=0;y<Ny;y++)
			{
				for(int m=0;m<7;m++) //�߽����� 
				{
					int a,b,c,d;
					a=x+1;
					b=x-1;
					c=y+1;
					d=y-1; 
					if(a==Nx)
					{
						a=1;
					}
					if(b<0)
					{
						b=Nx-2;//
					}
					if(c==Ny)
					{
						c=1;
					}
					if(d<0)
					{
						d=Ny-2;
					}
					delta_etas[m][x][y]=(etas[m][a][y]+etas[m][b][y]+etas[m][x][c]+etas[m][x][d]-4*etas[m][x][y])/(delta_x*delta_x);
				}
			}
		}

		
		for(int i=0;i<N;i++)
		{
			for(int x=0;x<Nx;x++)
			{
				for(int y=0;y<Ny;y++)
				{
					double sum=0.0;
					int n=0;
					for(int j=0;j<N;j++)
					{
						if(etas[j][x][y]!=0)
						{
	                        sum=sum+2*Mphi*(W*(etas[j][x][y]-etas[i][x][y])+a*a*(delta_etas[j][x][y]-delta_etas[i][x][y])/2-8*delta_E*sqrt(etas[i][x][y]*etas[j][x][y])/Pi);
							n=n+1;
		                    }
					//sum=sum+2*Mphi*(W*(etas[j][x][y]-etas[i][x][y])+a*a*(delta_etas[j][x][y]-delta_etas[i][x][y])/2-8*delta_E*sqrt(etas[i][x][y]*etas[j][x][y])/Pi)/25.;
					}
					//etas1[i][x][y]=-sum*delta_t+etas[i][x][y];
					etas1[i][x][y]=-sum*delta_t/n+etas[i][x][y];
				}
			}
		}
		
		
		for(int i=0;i<N;i++)
		{
			for(int x=0;x<Nx;x++)
			{
				for(int y=0;y<Ny;y++)
				{
					if(etas1[i][x][y]<=1.0e-5)
					{
						etas1[i][x][y]=0.0;
					}	
				}
			}
		}
		
		
		
		
		for(int x=0;x<Nx;x++)
		{
			for(int y=0;y<Ny;y++)
			{
				double sum=0.0;
				for(int i=0;i<7;i++)
				{
					sum=sum+etas1[i][x][y];
				}
				if ( sum!=0.0)
				{
					for(int i=0;i<N;i++)
					{
						etas[i][x][y]=etas1[i][x][y]/sum;
					}
				}
				
			}
		}
		double tttt,temp;
		tttt=time/TIME*100;
		temp=fmod(tttt,10);
		if(temp==0)
		{
			if(tttt==100)
			{
				cout<<".....������ɣ����VTK�ļ�.....";
			} 
			else
			{
				cout<<"�����У������"<<"\t"<<tttt<<".0%"<<endl; 
			}
			char fname[256];
			int a;
			a=tttt/10;
			double **PHI=new double *[Nx];
			for(int i=0;i<Nx;i++)
			{
				PHI[i]= new double[Ny];
			}
			for(int i=0;i<Nx;i++)
			{
				for(int j=0;j<Ny;j++)
				{
					double sum=0;
					for(int n=0;n<N;n++)
					{
						sum=sum+etas[n][i][j]*etas[n][i][j];
					}
					PHI[i][j]=sum;
				}
			}
			VTK(PHI,Nx,Ny,a,dx,dy);
		}
	}
//	fclose(in);//�����ļ� 
	return 0;
}
