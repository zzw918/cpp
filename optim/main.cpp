// 6个晶粒优化长大
#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <direct.h>
#include <windows.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#define PI 3.1415926
using namespace std;

// 多数情况下每一组的晶粒数是不会超过3的，只有在最终多个晶体交接处，才会出现晶粒数为4的情况。
// 注意：大于1的数字，必须要修正；小于10e-5的数字，必须要修正

// 创建晶粒结构体，key为晶粒标号，phi为对应的相场值。
struct Grain
{
    int key;
    double phi;
};

// 降序排序函数
bool LessSort(Grain a, Grain b)
{
    return (a.phi > b.phi);
};

const int N = 7;    // the total number of grains
const int Nx = 200; // the x-axis grid numbers
const int Ny = 120; // the y-axis grid numbers

int ti(int); // deal with x-axis periodic boundry
int tj(int); // deal with y-axis periodic boundry

// 输出函数
void output(vector<struct Grain> grid_b[][Ny], vector<struct Grain>::iterator it_b[][Ny], int N, int Nx, int Ny, int fileNum); // difine the output function
// 初始化相关
void init(Grain grain, vector<struct Grain> grid[][Ny], vector<struct Grain> grid_b[][Ny], int n, int i, int j);
// 获取特定序号(key)的phi值
double phi(vector<struct Grain> grid[][Ny], vector<struct Grain>::iterator it[][Ny], int key, int i, int j);

int main()
{
    cout << "Begin to calculate six grains with kim algorithm!" << endl;
    // grid为当前时刻网格vector，gird_b为下一时刻；it和it_b为对应的迭代器。
    vector<struct Grain> grid[Nx][Ny];
    vector<struct Grain> grid_b[Nx][Ny];

    vector<struct Grain>::iterator it[Nx][Ny];
    vector<struct Grain>::iterator it_b[Nx][Ny];

    // 演化方程运算中使用，即可相邻晶粒进行运算时使用
    vector<struct Grain>::iterator k[Nx][Ny];

    int round = 2;
    // 首先，对基体进行初始化
    Grain grain;
    grain = {6, 1.0}; // 创建标号为6的晶粒的phi值为1.0
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            // 所有格点都存在phi值为1的基体晶粒
            grid[i][j].push_back(grain);
            grid_b[i][j].push_back(grain);

            // 迭代器指向vector方便后续调用
            it[i][j] = grid[i][j].begin();
            it_b[i][j] = grid_b[i][j].begin();
        }
    }

    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            // 初始化必然涉及到每一个晶粒
            for (int n = 0; n < N; n++)
            {
                switch (n)
                {
                case 0:
                    if ((i <= Nx / 6.0 + round && i >= Nx / 6.0 - round) && (j <= Ny / 2.0 + round && j >= Ny / 2.0 - round))
                    {
                        init(grain, grid, grid_b, n, i, j);
                    }
                    break;

                case 1:
                    if ((i <= (Nx / 6.0) * 3 + round && i >= (Nx / 6.0) * 3 - round) && (j <= Ny / 2.0 + round && j >= Ny / 2.0 - round))
                    {
                        init(grain, grid, grid_b, n, i, j);
                    }
                    break;

                case 2:
                    if ((i <= (Nx / 6.0) * 5 + round && i >= (Nx / 6.0) * 5 - round) && (j <= Ny / 2.0 + round && j >= Ny / 2.0 - round))
                    {
                        init(grain, grid, grid_b, n, i, j);
                    }
                    break;

                case 3:
                    if ((i <= Nx / 3.0 + round && i >= Nx / 3.0 - round) && (j <= round || j >= Ny - round))
                    {
                        init(grain, grid, grid_b, n, i, j);
                    }
                    break;

                case 4:
                    if ((i <= (Nx / 3.0) * 2 + round && i >= (Nx / 3.0) * 2 - round) && (j <= round || j >= Ny - round))
                    {
                        init(grain, grid, grid_b, n, i, j);
                    }
                    break;

                case 5:
                    if ((i <= round && j <= round) || (i >= Nx - round && j <= round) || (i >= Nx - round && j >= Ny - round) || (i <= round && j >= Ny - round))
                    {
                        init(grain, grid, grid_b, n, i, j);
                    }
                    break;
                }
            }
        }
    }

    // 初始化相场输出
    output(grid_b, it_b, N, Nx, Ny, 0);

    // set the interval time and the whole time
    double deltaT = 0.01, // timeInterval
        allTime = 18.0;   // the whole time to grow

    double garma = 0.208, //  J/m2
        deltaE = 0.09,
           deltaX = 0.5e-6,
           Qb = 110e3,       //  j/mol
        R = 8.314,           //  j/(K*mol)
        T = 800.0,           //  T is Kelvin's temperature
        thigma = 7 * deltaX, // delta x is 0.5um，use the 'm'
        W = 4 * garma / thigma,
           a = (2 / PI) * pow(2 * thigma * garma, 0.5),
           M = (0.139 / T) * exp(-Qb / (R * T)) * PI * PI / (8 * thigma);

    double curTime = 0.0;
    int aid = 0;    // to help show output
    int number = 1; // filename number mark

    int files = 8; // file numbers
    int times = 0; // 计算每个格点计算演化方程的平均次数。
    while (curTime < allTime)
    {
        if ((aid % (int(allTime / deltaT) / files)) == 0)
        {
            cout << double(curTime / allTime) * 100 << "% has been calculated..." << endl;
        }

        aid++;

        // itKey表示当前晶粒的序号，通过序号寻找相应Phi值。
        int itKey = 0;
        // kKey表示相邻晶粒的序号，通过序号寻找相应Phi值。
        int kKey = 0;

        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                // 利用recodeLoc记录当前迭代器位置。
                int recordLoc = 0;
                int grainNum = 0;
                for (it[i][j] = grid[i][j].begin(), it_b[i][j] = grid_b[i][j].begin(); it[i][j] != grid[i][j].end(); it[i][j]++, it_b[i][j]++)
                {
                    double temp = 0.0;
                    double dif = 0.0;
                    itKey = it[i][j]->key;
                    times++;
                    for (k[i][j] = grid[i][j].begin(); k[i][j] != grid[i][j].end(); k[i][j]++)
                    {
                        double E = 0.0;
                        kKey = k[i][j]->key;
                        if (itKey == N - 1 && kKey != N - 1)
                        {
                            E = -0.09e6;
                        }
                        if (itKey != N - 1 && kKey == N - 1)
                        {
                            E = +0.09e6;
                        }
                        dif = ((phi(grid, it, kKey, ti(i + 1), j) + phi(grid, it, kKey, ti(i - 1), j) + phi(grid, it, kKey, i, tj(j + 1)) + phi(grid, it, kKey, i, tj(j - 1)) - 4 * phi(grid, it, kKey, i, j)) - (phi(grid, it, itKey, ti(i + 1), j) + phi(grid, it, itKey, ti(i - 1), j) + phi(grid, it, itKey, i, tj(j + 1)) + phi(grid, it, itKey, i, tj(j - 1)) - 4 * phi(grid, it, itKey, i, j))) / pow(deltaX, 2);
                        temp += (2 * M) * (W * (phi(grid, it, kKey, i, j) - phi(grid, it, itKey, i, j)) + 0.5 * pow(a, 2) * dif - 8 / PI * pow(phi(grid, it, itKey, i, j) * phi(grid, it, kKey, i, j), 0.5) * E);
                    }
                    // 因为这里的迭代器和vector都是全局的，所以迭代器被改变，下面进行修正
                    it[i][j] = grid[i][j].begin();   // 先指向最开始的位置
                    it[i][j] = it[i][j] + recordLoc; // 这样，it[i][j]的指向才会正确。
                    it_b[i][j]->phi = (-temp * deltaT / grid[i][j].size()) + it[i][j]->phi;

                    // phi值不会大于1，若大于1，则修正为1。
                    if (it_b[i][j]->phi > 1)
                    {
                        it_b[i][j]->phi = 1;
                    }

                    // phi值小于10e-5，则使之为0
                    if (it_b[i][j]->phi < 10e-5)
                    {
                        it_b[i][j]->phi = 0;
                    }

                    // remove phi值不会小于0，小于0的删去！
                    if ((it[i][j]->phi == 0) && ((-temp * deltaT / grid[i][j].size()) < 0))
                    {
                        // 这里必须是it_b[i][j]，而不是it_b[i][j]，指向不能混乱！ 因为之前设定了it_b[i][j] = grid_b[i][j].begin(); 而不是it[i][j]。
                        grid_b[i][j].erase(it_b[i][j]);
                        it_b[i][j]--;
                    }
                    recordLoc++; // recordLoc使得it[i][j]回到中断处
                    grainNum++;
                }

                // 先判断晶粒总数是否大于3(我们指定的最大晶粒数)，如果大于则排序，否则不用排序！
                if (grid_b[i][j].size() > 3)
                {
                    // 排序(从大到小)
                    sort(grid_b[i][j].begin(), grid_b[i][j].end(), LessSort);
                    // 删除前三大的结构体
                    int tem = 0;
                    for (it_b[i][j] = grid_b[i][j].begin(); it_b[i][j] != grid_b[i][j].end(); it_b[i][j]++)
                    {
                        if (tem > 2)
                        {
                            grid_b[i][j].erase(it_b[i][j]);
                            it_b[i][j]--;
                        }
                        tem++;
                    }
                    // 归一化处理
                    double sum = 0.0;
                    for (it_b[i][j] = grid_b[i][j].begin(); it_b[i][j] != grid_b[i][j].end(); it_b[i][j]++)
                    {
                        sum += it_b[i][j]->phi;
                    }

                    for (it_b[i][j] = grid_b[i][j].begin(); it_b[i][j] != grid_b[i][j].end(); it_b[i][j]++)
                    {
                        double temp = it_b[i][j]->phi;
                        it_b[i][j]->phi = temp / sum;
                    }
                }

                // 添加符合要求的晶粒进入网格结构体中。 其中own存储已有额key， oth存储不存在的key，两者都是vector。
                vector<int> own;
                vector<int>::iterator itown;
                for (it_b[i][j] = grid_b[i][j].begin(); it_b[i][j] != grid_b[i][j].end(); it_b[i][j]++)
                {
                    own.push_back(it_b[i][j]->key);
                }

                vector<int> oth;
                vector<int>::iterator itoth;
                //N - 1，是因为基体晶粒只可能消失，而不可能扩张，所以不加
                for (int i = 0; i < N - 1; i++)
                {
                    bool flag = false;
                    for (itown = own.begin(); itown != own.end(); itown++)
                    { // 2, 3, 5
                        if (i == *itown)
                        {
                            flag = true;
                            break;
                        }
                    }
                    if (flag == false)
                    {
                        oth.push_back(i);
                    }
                }

                // oth是需要在其他的网格中(最近的四个格点)寻找的晶粒，注意，这里的寻找，应该使用的是当前时刻grid[i][j]，而不是下一时刻grid_b[i][j]
                // 因为我们将之将入下一时刻，是因为当前时刻可能会受到影响！所以添加之后，总数可能就会大于3了！
                for (itoth = oth.begin(); itoth != oth.end(); itoth++)
                {
                    bool flag = false;
                    // 右方格点
                    for (it[ti(i + 1)][j] = grid[ti(i + 1)][j].begin(); it[ti(i + 1)][j] != grid[ti(i + 1)][j].end(); it[ti(i + 1)][j]++)
                    {
                        if (it[ti(i + 1)][j]->key == *itoth && it[ti(i + 1)][j]->phi != 0)
                        {
                            grain = {*itoth, 0.0};
                            // 这里的添加使用的时grid_b[i][j]，而不是grid[i][j]
                            grid_b[i][j].push_back(grain);
                            // 添加之后，我们就可以跳过这个晶粒，去看其他标号的晶粒是否满足条件，即flag = true;以下同理，不再赘述。
                            flag = true;
                            break;
                        }
                    }
                    if (flag == true)
                    {
                        continue;
                    }

                    // 左方格点
                    for (it[ti(i - 1)][j] = grid[ti(i - 1)][j].begin(); it[ti(i - 1)][j] != grid[ti(i - 1)][j].end(); it[ti(i - 1)][j]++)
                    {
                        if (it[ti(i - 1)][j]->key == *itoth && it[ti(i - 1)][j]->phi != 0)
                        {
                            grain = {*itoth, 0.0};
                            grid_b[i][j].push_back(grain);
                            flag = true;
                            break;
                        }
                    }
                    if (flag == true)
                    {
                        continue;
                    }

                    // 上方格点
                    for (it[i][tj(j + 1)] = grid[i][tj(j + 1)].begin(); it[i][tj(j + 1)] != grid[i][tj(j + 1)].end(); it[i][tj(j + 1)]++)
                    {
                        if (it[i][tj(j + 1)]->key == *itoth && it[i][tj(j + 1)]->phi != 0)
                        {
                            grain = {*itoth, 0.0};
                            grid_b[i][j].push_back(grain);
                            flag = true;
                            break;
                        }
                    }
                    if (flag == true)
                    {
                        continue;
                    }

                    // 下方格点
                    for (it[i][tj(j - 1)] = grid[i][tj(j - 1)].begin(); it[i][tj(j - 1)] != grid[i][tj(j - 1)].end(); it[i][tj(j - 1)]++)
                    {
                        if (it[i][tj(j - 1)]->key == *itoth && it[i][tj(j - 1)]->phi != 0)
                        {
                            grain = {*itoth, 0.0};
                            grid_b[i][j].push_back(grain);
                            break;
                        }
                    }
                }
            }
            // 深克隆赋值不能在这里，否则，与当前格点(i, j)相邻的格点计算会出现错误！
        }

        // 当所有的格点都进行运算之后，进行深克隆赋值，保证后续时间步的晶粒一一对应。
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                grid[i][j].clear(); // 清空grid[i][j]，否则无法对应
                for (it_b[i][j] = grid_b[i][j].begin(); it_b[i][j] != grid_b[i][j].end(); it_b[i][j]++)
                {
                    grain = {it_b[i][j]->key, it_b[i][j]->phi};
                    grid[i][j].push_back(grain); // 完成赋值
                }
                it[i][j] = grid[i][j].begin(); // 绑定iterator
            }
        }

        // time added
        curTime += deltaT;

        if (aid % ((int)(allTime / deltaT) / files) == 0)
        {
            output(grid_b, it_b, N, Nx, Ny, number++);
        }
    }
    cout << "calculate end!" << endl
         << endl;
    system("pause");
    return 0;
}

// 初始化相关函数
void init(Grain grain, vector<struct Grain> grid[][Ny], vector<struct Grain> grid_b[][Ny], int n, int i, int j)
{
    grain = {n, 1.0};
    grid[i][j].push_back(grain);
    grid_b[i][j].push_back(grain);

    // 删除基体晶粒结构体
    grid[i][j].erase(grid[i][j].begin());
    grid_b[i][j].erase(grid_b[i][j].begin());
}

// 处理演化方程相关的计算
double phi(vector<struct Grain> grid[][Ny], vector<struct Grain>::iterator it[][Ny], int key, int i, int j)
{
    // 如果可以找到该晶格内的相应晶粒，就返回相场值，否则，就返回0
    for (it[i][j] = grid[i][j].begin(); it[i][j] != grid[i][j].end(); it[i][j]++)
    {
        if (it[i][j]->key == key)
        {
            return it[i][j]->phi;
        }
    }
    return 0.0;
}

// define the output function
void output(vector<struct Grain> grid_b[][Ny], vector<struct Grain>::iterator it_b[][Ny], int N, int Nx, int Ny, int fileNum)
{
    double phi_o[Nx][Ny];
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            double sum = 0.0;
            for (it_b[i][j] = grid_b[i][j].begin(); it_b[i][j] != grid_b[i][j].end(); it_b[i][j]++)
            {
                sum = sum + pow(it_b[i][j]->phi, 2);
            }
            phi_o[i][j] = sum;
        }
    }

    // output the vtk files
    fstream outfile;
    char filename[20] = "output";
    char num[20];
    char foldername[20] = "./outputfolder/";
    // 不能使用itoa函数，因为其不是标准函数，c++11规范不支持，所以这里使用sprintf
    sprintf(num, "%d", fileNum);
    strcat(strcat(filename, num), ".vtk");
    outfile.open(strcat(foldername, filename), ios::out);

    // 异常处理
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

// ti和tj都是为了处理周期边界条件所设定的。
int ti(int i)
{
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
