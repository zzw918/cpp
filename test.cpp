#include <iostream>
#include <vector>
#include <windows.h>
using namespace std;
int main() {
    // 1、使用vector
    // vector<int> test;
    // test.push_back(1);
    // test.push_back(2);
    // cout << test[0] << endl;

    // 2、使用共同体union
    // union icf {
    //     int i;
    //     double j;
    // };
    // icf arrDate[10];
    // arrDate[0].i = 1;
    // arrDate[0].j = 3.14;
    // arrDate[1].i = 2;
    // arrDate[1].j = 6.28;
    
    // 3、数组指针
    // vector <int> test;
    // test.push_back(55);
    // test.push_back(44);
    // test.push_back(22);
    // vector <int> * p = &test;
    // cout << (*p)[0] << endl;

    // icf phi[10][10][10];
    // phi[0][0][0].i = 5;
    // phi[0][0][0].j = 3.14;
    // cout << phi[0][0][0].i << endl << phi[0][0][0].j << endl;
    // phi[10][10] = icf;    
    // for (int i = 0; i < 10; i++) {
    //     for (int j = 0; j < 10; j++) {
    //         phi[i][j] = &arrDate;
    //     }
    // }

    // 建立一个格点的共同体，每个格点上存储着有phi值的晶粒标号以及其对应的phi值。
    struct grains
    {
        int i;
        double j;
    }; 
    // union grains {
    //     int i; // i是晶粒的标号
    //     double j; // j是i晶粒的phi值
    // };
    // grains point1[10];
    // point1[0].i = 5;
    // point1[0].j = 0.2; // 在point1格点上的第一个晶粒
    // point1[1].i = 8;
    // point1[1].j = 0.5; // 在point1格点上的第二个晶粒
    // point1[2].i = 20;
    // point1[2].j = 0.3; // 在point1格点上的第三个晶粒

    // cout << point1[0].i << endl << point1[1].j << endl << point1[2].i << endl; // 5 0.5 20

    grains* phi[20][20]; // 建立一个20 * 20的格点数组，每一个格点都指向一个grains结构体
    // 使得每一个格点，指向一个具体的结构体
    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            grains point[10]; 
            point[0].i = i;
            point[0].j = 0.1;
            phi[i][j] = &point[0];
        }
    }
    cout << (*phi[1][5]).i << endl; // 19

    // grains* pPhi;
    // grains point[10];
    // point[0].i = 2;
    // point[0].j = 0.1;
    // point[1].i = 4;
    // point[1].j = 0.2;
    // pPhi = &point[0];
    // cout << (*(pPhi+1)).i << endl; // 2
    system("pause");
}