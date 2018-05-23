// 6进行动态再结晶
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

// ATTENTION: 再结晶阶段，储能(驱动力)E是否始终为0，还是将最开始长大后的六边形晶粒类似于之前的基体；如果为0，那么就相当于晶粒粗话，大的增大，小的减小，这样，刚生长出来的晶体
// 存活不就就会消失；如果是把六边形晶粒当做是基体，储能(驱动力)参与计算，那么E值应该不是静止不动的。
// E应该是变化的，否则，晶粒长大一部分之后，就会消失。

// ATTENTION: 晶粒再结晶过程中，达到临界位错密度需要较长的时间，所以，我们可以把再结晶的时间继续放大。
// FIXME: 再结晶过程中，为什么还会有大于0的值(1.1)呢？能否强行修正，在output函数中！我觉得是可以的！

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

// 本程序为晶粒再结晶程序。在晶粒长大阶段，N值不会发生变化；但在再结晶时，会在晶界处生成新的晶核，所以N值会发生变化。主要用于序号。
int N = 7;    
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
// 获取(i, j)中Key晶粒的位错密度。
double getRho(int i, int j, vector<struct Grain> grid[][Ny], vector<struct Grain>::iterator itNew[][Ny], double rho_b[][Ny], int Key);

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

    // 获取rho值时使用。
    vector<struct Grain>::iterator itNew[Nx][Ny];

    int round = 5;
    // 首先，对基体进行初始化
    Grain grain;
    grain = {N - 1, 1.0}; // 创建标号为6的晶粒的phi值为1.0
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
    double deltaT = 0.01; // timeInterval

    double garma = 0.208, //  J/m2
           allTime = 25.0,   // the whole time to grow    
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
    // 结束之后，grid和grid_b的值是保持的，可以继续在下面的再结晶过程中使用。
    cout << "calculate end, grain grown, begin to dynamic recrystallization!" << endl << endl;

    cout << "DRX..." << endl;
    // 以上，晶粒长大结束，下面进行动态再结晶过程的计算。
    // 动态再结晶是基于长大的过程进行计算的
    // 这里需要使用到一些新的变量，如rho等, A1, A2 
    double rho0 = 10e9, // 初始位错密度
           //sigama_c = 4.0e7, // 40Mpa
           //rhoc = pow((sigama_c / (A * U * B)), 2), 
           rhoc = 5.51e13, // 临界位错密度
           k1 = 4e8, // 许婷
           A = 0.5, // 即aerfa
           A1 = 2.0e44,
           A2 = 7.6, // A1和A2用于计算sigma_s
           U = 4.21e10, // 即μ
           B = 2.56e-10, // 即伯氏矢量b
           epsilon = 2e-3, //即埃普西隆，小e 许婷
           Qa = 2.75e5, // 激活能
           c = 5.0e25, // d和c是常量
           d = 1.0,
           sigma_s = pow(A1 * epsilon * exp(Qa / (R * T)), 1/A2), // 许婷
           k2 = A * U * B * k1 / sigma_s,
           delta_n = c * pow(epsilon, d) * exp(-Qa / (R * T)); /// 许婷
           // 既然 ngb 和 step都是变化的，所以放在这里是不合适的，应该放在需要用到的位置。
    
    // 驱动力相关参数
    double tao = 0.5 * U * B  * B;
    // E = tao * (rhoi - rhoj); // 这里是驱动力E的计算方式

    // cout << "sigama_c" << sigama_c << endl
    //      << "sigma_s" << sigma_s << endl
    //      << "delta_n" << delta_n << endl
    //      << "k2" << k2 << endl;
    //     //  每个格点上都存储位错密度，值为rho0，double类型，因为此时还没有加应变，所以认为各处的位错密度相等。

    static double rho[Nx][Ny]; // 上一时刻的位错密度
    static double rho_b[Nx][Ny]; // 当前时刻的位错密度

    // 设置参考点，减少后续循环的计算！
    double refPoint = rho0;
    double refPointB = rho0;

    // 设置初始位错密度
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            rho[i][j] = rho0;
            rho_b[i][j] = rho0;
        }
    }

    // 用于存储境界上的晶粒格点位置
    // 注意：这里要开辟的空间比较大，所以，我们需要使用static全局变量。
    static int pos[Nx * Ny][2];

    double curTimeDRX = 0.0; // 时间从新开始
    int aidDRX = 0;    // to help show output
    int numberDRX = 100; // 标号从晶粒长大末尾更后面一些的标号开始，先防止覆盖，优化之后再做
    double allTimeDRX = 40.0;

    int maxNum = 5; // 同一网格内phi值不为0的晶体个数，用于排序算法

    int filesDRX = 8; // file numbers

    int stepFlag = 0;

    int ngb = 0; // 这里ngb记录的满足条件（rho达到rhoc并且在境界上）的网格个数    
    // step是根据ngb的不同，而发生变化的，即要求的时间步
    // double step = pow(delta_n * deltaT * ngb * deltaX * deltaX / thigma, -1); // 许婷
    // 开始ngb一定为0，所以，这里计算会出错，step为一个不可计算值
    double step = 0;
    

    // 再结晶演化计算
    while (curTimeDRX < allTimeDRX)
    {
        if ((aidDRX % (int(allTimeDRX / deltaT) / filesDRX)) == 0)
        {
            cout << double(curTimeDRX / allTimeDRX) * 100 << " % has been calculated..." << endl;
        }

        //循环次数，输出文件使用
        aidDRX++;

        // 每个事件步开始，都需要施加应力，进而位错密度发生变化 --- 这个变化一定是始终增加的！所以一旦第一次发生了再结晶，后面每一步都会再结晶。
        // 在这个step的最后，需要将当前rho_b赋值给上一时刻位错密度数组rho，以保证循环使用
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                rho_b[i][j] = (k1 * pow(rho[i][j], 0.5) - k2 * rho[i][j]) * (epsilon * deltaT) + rho[i][j];
            }
        }
        // ATTENTION: 相对参考点同时计算！
        refPointB = (k1 * pow(refPoint, 0.5) - k2 * refPoint) * (epsilon * deltaT) + refPoint;

        // 遍历查找有多少处在边界，并且符合位错密度条件的
        // 1、因为边界的比较少，所以我们先判断是否在边界，后判断是否达到了临界密度
        // 2、但是如果先判断边界，则需要迭代器遍历每一个格点，这个是非常耗时的；而比较位错密度却相对容易一些，计算复杂度小。
        // 综上，选择第二种方案。

        // 如果达到了下面两个要求，就可以生成新的晶粒。
        // 1、已经隔了特定的时间步(这里参考英文文献中时间步的说法，比如每隔10步生成一个新的晶粒)
        // 2、ngb的值是大于0的，这样，我们才能在晶界上任意选择一个区域生成一个晶核
        // 下面stepFlag = step + 1时生成新的晶粒(或者是stepFlag > step，结果相同)

        // 在refPointB > rhoc之后，就会不断产生新的晶核； 而在此之前，不会产生新的晶核，那么之前的六边形不会变，所以，不用进行演化方程的计算。
        if (refPointB > rhoc) {
            //step次数，产生新的晶粒使用，只需在refPintB > rhoc之后计算，之前的统计都是没有意义的。
            stepFlag++;
            if (stepFlag > step)
            {
                // 每次进入可以生成晶核的节点，ngb都是要计算的。
                ngb = 0; //ATTENTION: ngb需要置0，否则会进入累加循环地狱
                for (int i = 0; i < Nx; i++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        if (rho_b[i][j] > rhoc)
                        {
                            double sum = 0.0;
                            for (it[i][j] = grid[i][j].begin(); it[i][j] != grid[i][j].end(); it[i][j]++)
                            {
                                sum += pow(it[i][j]->phi, 2);
                            }
                            if (sum < 0.6)
                            {
                                // 每隔step时间步，也是不需要清空这个数组的，因为我们只要能准确记录ngb最终大小以及对应的坐标即可
                                pos[ngb][0] = i;
                                pos[ngb][1] = j;
                                // 如果存在，则ngb++，如果从0开始，那么第一个符合之后，ngb的值就是1了。
                                ngb++;
                            }
                        }
                    }
                }
    
                if (ngb > 0) {
                    // 生成新的晶粒
                    srand(time(0));
                    N = N + 1;            // 即此时晶粒个数增加，所以将程序开始出的N从const修改为一般的N值。
                    grain = {N - 1, 1.0}; // 创建新的晶粒结构体，因为标号从0开始，所以，这里key是N - 1，因为是添加晶核，phi为1.0
                    cout << "new grain" << endl;
                    // 因为在晶界处任意放置，所以，我们先随机选取一个点。
                    int index = rand() % ngb; // 这里得到的序号一定是随即的，并且是符合要求的
                    for (int i = 0; i < Nx; i++)
                    {
                        for (int j = 0; j < Ny; j++)
                        {
                            int dif = pow((i - pos[index][0]), 2) + pow((j - pos[index][1]), 2) - pow(round, 2);
                            if (dif < 0 || dif == 0)
                            {
                                // 在某一格点处，一定是其他晶粒phi值均不存在，所以首先清空
                                grid[i][j].clear();
                                grid_b[i][j].clear();
    
                                // 然后将当前晶粒的grain储存在grid[i][j]中。
                                grid[i][j].push_back(grain);
                                grid_b[i][j].push_back(grain);
    
                                // 接着把其rho设置为rho0
                                rho[i][j] = rho0;
                                rho_b[i][j] = rho0;
                            }
                        }
                    }
                    // step是根据ngb的不同，而发生变化的，即要求的时间步。如果ngb > 0，这里正常计算；
                    step = pow(delta_n * deltaT * ngb * deltaX * deltaX / thigma, -1);
                } 
                else 
                {
                    // 如果ngb == 0，则step计算就会出错，即1 / 0，所以，我们直接将step设置为0，下一步继续可以进入计算！
                    step = 0;
                }
    
                // 这里不管是否满足条件，都需要成为0，进行下一次运算。
                stepFlag = 0;
            }
    
    
            // 下面这一部分和晶粒长大部分使用的程序大致相同，只是在驱动力的计算上有些区别，其他一致。 
    
            // itKey表示当前晶粒的序号，通过序号寻找相应Phi值。
            int itKey = 0;
            // kKey表示相邻晶粒的序号，通过序号寻找相应Phi值。
            int kKey = 0;

            // 这里的rhoi和rhoj是在计算E(驱动力)需要的不同位置的位错密度。
            double rhoi = 0;
            double rhoj = 0;
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    // 利用recodeLoc记录当前迭代器位置。
                    int recordLoc = 0;
                    int grainNum = 0;
                    // int gridSize = 0;
                    for (it[i][j] = grid[i][j].begin(), it_b[i][j] = grid_b[i][j].begin(); it[i][j] != grid[i][j].end(); it[i][j]++, it_b[i][j]++)
                    {
                        // FIXME: 再结晶时，已经不存在基体了，还需要使用E吗？或者说，我们需要让最初的几个晶粒类似于之前的基体？ 这个需要实验。
                        // 暂时是需要使用的!
                        double temp = 0.0;
                        double dif = 0.0;
                        itKey = it[i][j]->key;
                        // m为横坐标，n为纵坐标；注意m、n的范围，这是为了减小复杂度计算。
                        // gridSize = grid[i][j].size();

                        // if (gridSize > 1) {
                            rhoi = getRho(i, j, grid, itNew, rho_b, itKey);
                            // cout << " rhoi: " << rhoi;
                        // }
                        
                        for (k[i][j] = grid[i][j].begin(); k[i][j] != grid[i][j].end(); k[i][j]++)
                        {
                            double E = 0.0;
                            kKey = k[i][j]->key;
                            // grid[i][j]格点晶粒个数只有两种情况：等于1或者大于1；如果等于1，那么E就是0，否则，就要进行计算！
                            // if (gridSize > 1)
                            // {
                                rhoj = getRho(i, j, grid, itNew, rho_b, kKey);
                                // cout << " rhoj: " << rhoj;
                                E = tao * (rhoj - rhoi);
                                // cout << " E: " << E;
                            // }
                            dif = ((phi(grid, it, kKey, ti(i + 1), j) + phi(grid, it, kKey, ti(i - 1), j) + phi(grid, it, kKey, i, tj(j + 1)) + phi(grid, it, kKey, i, tj(j - 1)) - 4 * phi(grid, it, kKey, i, j)) - (phi(grid, it, itKey, ti(i + 1), j) + phi(grid, it, itKey, ti(i - 1), j) + phi(grid, it, itKey, i, tj(j + 1)) + phi(grid, it, itKey, i, tj(j - 1)) - 4 * phi(grid, it, itKey, i, j))) / pow(deltaX, 2);
                            temp += (2 * M) * (W * (phi(grid, it, kKey, i, j) - phi(grid, it, itKey, i, j)) + 0.5 * pow(a, 2) * dif - 8 / PI * pow(phi(grid, it, itKey, i, j) * phi(grid, it, kKey, i, j), 0.5) * E);
                        }
                        // 因为这里的迭代器和vector都是全局的，所以迭代器被改变，下面进行修正
                        // 补充：在函数传递时，如果在函数中新建迭代器，则不需要修正，后续优化。
                        it[i][j] = grid[i][j].begin();   
                        it[i][j] = it[i][j] + recordLoc; 
                        
                        // 获得新时刻的phi值 
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
                    if (grid_b[i][j].size() > maxNum)
                    {
                        // 排序(从大到小)
                        sort(grid_b[i][j].begin(), grid_b[i][j].end(), LessSort);
                        // 删除前三大的结构体
                        int tem = 0;
                        for (it_b[i][j] = grid_b[i][j].begin(); it_b[i][j] != grid_b[i][j].end(); it_b[i][j]++)
                        {
                            if (tem > maxNum - 1)
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
    
                    // 添加符合要求的晶粒进入网格结构体中。 其中own存储已有的key， oth存储不存在的key，两者都是vector。
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
        }

        // 将rho_b赋值给rho，以便下一次循环中使用。
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                rho[i][j] = rho_b[i][j];
            }
        }

        // 同样地，将参考坐标点也赋值给上一时刻。
        refPoint = refPointB;

        // time added
        curTimeDRX += deltaT;

        if (aidDRX % ((int)(allTimeDRX / deltaT) / filesDRX) == 0)
        {
            output(grid_b, it_b, N, Nx, Ny, numberDRX++);
        }
    }
    cout << "Finished! " << endl;
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
            // FIXME: 这里可以强行修正吗？
            if (sum > 1) {
                sum = 1;
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


// FIXME: 注意超出边界的问题! 会出现错误！
// 但是，如果不使用这种方法，效率就会非常慢，所以，这里处理好边界问题是当务之急！
double getRho(int i, int j, vector<struct Grain> grid[][Ny], vector<struct Grain>::iterator itNew[][Ny], double rho_b[][Ny], int Key)
{
    // 这里取 Nx/10 以及 Ny/10 的范围是足够的!
    // for (int m = i - Nx / 10; m < i + Nx / 10; m++)
    // {
    //     for (int n = j - Ny / 10; n < j + Ny / 10; n++)
    //     {
    //         for (itNew[m][n] = grid[m][n].begin(); itNew[m][n] != grid[m][n].end(); itNew[m][n]++)
    //         {
    //             if (itNew[m][n]->key == Key)
    //             {
    //                 if ( && itNew[m][n]->phi > 1.0) {
    //                     cout << " rho_b: " << rho_b[m][n];
    //                     return rho_b[m][n];
    //                 }
    //             }
    //         }
    //     }
    // }
    for (int m = 0; m < Nx; m++)
    {
        for (int n = 0; n < Ny; n++)
        {
            // 第一种方法：
            for (itNew[m][n] = grid[m][n].begin(); itNew[m][n] != grid[m][n].end(); itNew[m][n]++)
            {
                if (itNew[m][n]->key == Key && itNew[m][n]->phi == 1.0)
                {
                    return rho_b[m][n];
                }
            }
            // 第二种方法：
            // grid[m][n]在上一轮的添加过程中，
            // if (grid[m][n][0].key == Key && grid[m][n][0].phi == 1.0)
            // {
            //     return rho_b[m][n];
            // }
        }
    }
}
