#include <iostream>
#include <windows.h>
#include <ctime>
#include <cmath>
#include <vector>
using namespace std;
struct GrainRho
{
    int key;
    double rho;
};

int main() {
    // 创建一个存储晶粒位错密度的结构体vector，结构体中存储key以及rho值：
    vector<struct GrainRho> rhovec;
    vector<struct GrainRho> rhovec_b;

    // 创建相应的迭代器
    vector<struct GrainRho>::iterator itrho;
    vector<struct GrainRho>::iterator itrho_b;

    itrho = rhovec.begin();
    itrho_b = rhovec.begin();



    system("pause");
}
