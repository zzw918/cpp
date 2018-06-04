// 随机生成晶粒并长大
#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <direct.h>
#include <windows.h>
#include <vector>
#include <ctime>
#include <algorithm>
#include <stdlib.h>
#define PI 3.1415926
using namespace std;

// 这个程序中，我们放了11个晶粒，包括一个基体晶粒，然后在长大的过程中，大的晶粒不断长大，而小的晶粒不断变小甚至消失，这就是晶粒粗化的过程
// 并且最终成功由11个晶粒变成了6个晶粒（理论上晶粒还会继续粗化，小的晶粒继续变小、消失），实现了晶粒粗化

// 注意：大于1的数字，必须要修正；小于10e-5的数字，必须要修正

// 创建晶粒结构体，key为晶粒标号，phi为对应的相场值。
struct Grain
{
    int key;
    double phi;
};

// 随机数结构体，存储符合结晶条件的位置的晶粒
struct Position
{
    int randX; // 随机x坐标 
    int randY; // 随机y坐标
    int key;
};


// 降序排序函数
bool LessSort(Grain a, Grain b)
{
    return (a.phi > b.phi);
};

// 随机晶粒，所以直接在一个方形区域内生成即可。25个晶粒，其中最后一个应该是基体，其他的我们认为是正常生长的。
const int N = 31;    // the total number of grains
const int Nx = 300; // the x-axis grid numbers
const int Ny = 200; // the y-axis grid numbers

// 用于记录时间
char tmp[64];

int ti(int); // deal with x-axis periodic boundry
int tj(int); // deal with y-axis periodic boundry

// 输出函数
void output(double phi_o[][Ny], vector<struct Grain> grid_b[][Ny], vector<struct Grain>::iterator it_b[][Ny], int N, int Nx, int Ny, int fileNum); // difine the output function
// 获取特定序号(key)的phi值
double phi(vector<struct Grain> grid[][Ny], vector<struct Grain>::iterator it[][Ny], int key, int i, int j);
char* getTime();

int main()
{
    // 创建文件，保存重要信息
    fstream outfile;
    char DFILE[30] = "./outputfolder/log.txt";
    outfile.open(DFILE, ios::out|ios::app);
    
    cout << "Begin to calculate random " << N - 1 << "-grains with kim algorithm!" << endl;
    cout << "Nx: " << Nx << ", Ny: " << Ny << "." << endl;
    
    // 通过outfile保存信息
    outfile << getTime() << endl;
    outfile << "Begin to calculate random " << N - 1 << "-grains with kim algorithm!" << endl;
    outfile << "Nx: " << Nx << ", Ny: " << Ny << "." << endl;

    // grid为当前时刻网格vector，gird_b为下一时刻；it和it_b为对应的迭代器。
    static vector<struct Grain> grid[Nx][Ny];
    static vector<struct Grain> grid_b[Nx][Ny];

    static vector<struct Grain>::iterator it[Nx][Ny];
    static vector<struct Grain>::iterator it_b[Nx][Ny];

    // 演化方程运算中使用，即可相邻晶粒进行运算时使用
    static vector<struct Grain>::iterator k[Nx][Ny];

    // 输出二维数组
    static double phi_o[Nx][Ny];

    // 首先，对基体进行初始化
    Grain grain;
    grain = {N - 1, 1.0}; // 创建标号为N-1的晶粒的phi值为1.0
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

    // randX为随机晶粒的横坐标
    // randY为随机晶粒的纵坐标
    // round为初始化晶粒的半径
    int randX = 0;
    int randY = 0;
    int round = 7;
    srand(time(0)); //生成随机数种子

    // 创建一维vector，用来存储N - 1对(所有的)随机坐标，并创建对应的迭代器
    static vector<struct Position> allPos;
    static vector<struct Position>::iterator itall;

    // 创建结构体实例，用于存储randX和randY
    Position position;    

    // 只需要初始化N - 1个随机晶粒，第N个晶粒是基体，已经初始化。
    for (int n = 0; n < N - 1; n++)
    {
        randX = rand() % Nx;
        randY = rand() % Ny;
        // 得到randX和randY之后，并不是直接就放在allPos里，因为可能多次得到的有重叠的情况，这样情况下，我们是不要的。
        bool flag = false;
        for (itall = allPos.begin(); itall != allPos.end(); itall++)
        {
            int dif = pow(pow((randX - itall->randX), 2) + pow((randY - itall->randY), 2), 0.5); 
            // dif得到的是当前点与已有点的直线距离，如果大于两倍半径，就成立；否则，重新获得
            if (dif < 2 * (round + 5)) {
                flag = true;
                break;
            } 
        }
        if (flag == false) {
            position = {randX, randY, n};
            allPos.push_back(position);
        } else {
            n--;
        }
    }

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int roundNum = 0;
            int realRound = round;
            for (itall = allPos.begin(); itall != allPos.end(); itall++)
            {
                // cout << "vector: " << itall->randX << " " << itall->randY << endl;
                // 一点是否在圆内，通过(x - x0)2 + (y - y0)2 <= r2来判断，即如果满足条件，说明(x, y)是在园内的，否则不再。
                roundNum++;
                if (roundNum % 2 == 0) {
                    realRound += 3;
                } else {
                    realRound -= 3;
                }
                int dif = pow(i - itall->randX, 2) + pow(j - itall->randY, 2) - pow(realRound, 2);
                if (dif < 0 || dif == 0) {
                    grain = {itall->key, 1.0};
                    // 添加这个晶粒到格点
                    grid[i][j].push_back(grain);
                    grid_b[i][j].push_back(grain);
                    // 删除基体晶粒结构体
                    grid[i][j].erase(grid[i][j].begin());
                    grid_b[i][j].erase(grid_b[i][j].begin());

                    // 如果满足条件并添加晶粒，那么其他的点就不用考虑了，直接break
                    break;
                }
            }
        }
    }

    // set the interval time and the whole time
    double deltaT = 0.01; // timeInterval

    double garma = 0.208, //  J/m2
        deltaE = 0.09,
        deltaX = 0.5e-6, // 单位是m
        Qb = 110e3,       //  j/mol
        R = 8.314,           //  j/(K*mol)
        T = 400.0,           //  T is Kelvin's temperature
        thigma = 7 * deltaX, // delta x is 0.5um，use the 'm'
        W = 4 * garma / thigma,
        a = (2 / PI) * pow(2 * thigma * garma, 0.5),
        M = (0.139 / T) * exp(-Qb / (R * T)) * PI * PI / (8 * thigma);

    double area = (Nx * deltaX) * (Ny * deltaX); // area为总面积，deltaX即每个格点的长度。

    double curTime = 0.0;
    double timeFileInterval = 5;
    double timeFlag = 0;
    int number = 1; // filename number mark


    cout << "temperature T is: " << T << endl;
    cout << "deltaT is " << deltaT << ", all files will be created every " << timeFileInterval / deltaT << " steps." << endl;

    // 通过outfile保存信息
    outfile << "temperature T is: " << T << endl;
    outfile << "deltaT is " << deltaT << ", all files will be created every " << timeFileInterval / deltaT << " steps." << endl << endl;
    outfile << "equivalent diameter: " << endl;

    // 关闭outfile，后续需要时再打开写入。
    outfile.close();

    cout << endl << "Initiated!" << endl;
    cout << "Please wait to computed..." << endl;
    // 相场初始化输出，但是这里需要先将前面的phi值设置为0
    output(phi_o, grid_b, it_b, N, Nx, Ny, 0);
    // system("pause");

    int files = 20; // file numbers
    int times = 0; // 计算每个格点计算演化方程的平均次数。

    while (true)
    {
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
                }
         
                // 归一化处理 - 此处理必须要做
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
        timeFlag += deltaT;

        if (timeFlag > timeFileInterval)
        {
            cout << number << " files have been created, current time is " << int(curTime) << endl;
            output(phi_o, grid_b, it_b, N, Nx, Ny, number++);
            timeFlag = 0;

            vector<int> grains;
            vector<int>::iterator itgrains;
            
            // 在这里我们可以继续获得晶粒的平均晶粒尺寸。即需要获得计算域总面积，再获得当前晶粒的总个数，最后通过公式即可计算得到所有晶粒的平均晶粒尺寸。
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    // 当前网格的晶粒数为1，表明此网格在晶粒内部，这是我们需要的。
                    if (grid[i][j].size() == 1) {
                        int key = grid[i][j][0].key;
                        bool flag = false;
                        for (itgrains = grains.begin(); itgrains != grains.end(); itgrains++)
                        {   
                            if (*itgrains == key) {
                                flag = true;
                                break;
                            }
                        }
                        if (flag == false) {
                            grains.push_back(key);            
                        }
                    }
                }
            }

            outfile.open(DFILE, ios::out|ios::app);
            outfile << getTime() << endl;
            // 这一部分内容可以建立一个log系统，然后将之进行记录，方便后续排查问题
            cout << "The number of grains: " << grains.size() << endl;
                
            double singleArea = area/grains.size(); // 获得每一个晶粒的面积大小
            double D = pow((singleArea / PI), 0.5) * 2; // 获得当量直径
            outfile << "number: " << number << ", curTime: " << int(curTime) << ", size: " << grains.size() << ", D: " << D << endl << endl;
            outfile.close();
        }
    }
    cout << "calculate end!" << endl
         << endl;
    system("pause");
    return 0;
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
void output(double phi_o[][Ny], vector<struct Grain> grid_b[][Ny], vector<struct Grain>::iterator it_b[][Ny], int N, int Nx, int Ny, int fileNum)
{
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            double sum = 0.0;
            for (it_b[i][j] = grid_b[i][j].begin(); it_b[i][j] != grid_b[i][j].end(); it_b[i][j]++)
            {
                sum = sum + pow(it_b[i][j]->phi, 2);
            }
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

char* getTime()
{
    time_t t = time(0);
    strftime(tmp, sizeof(tmp), "%Y/%m/%d %X", localtime(&t));
    return tmp;
}
