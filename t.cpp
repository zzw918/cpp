#include <iostream>
#include <windows.h>
#include <vector>
#include <algorithm>
using namespace std;
// 结构体声明
struct Grain {
    int key; 
    double phi;
};

const int N = 7;
const int Nx = 200;
const int Ny = 200;

// 排序函数
bool LessSort(Grain a, Grain b) {
    return (a.phi > b.phi);
};

int main() {
    // 结构体vector
    vector<struct Grain> grid[Nx][Ny];
    vector<struct Grain> grid_b[Nx][Ny];
    // 结构体vector迭代器
    vector<struct Grain>::iterator it[Nx][Ny];
    vector<struct Grain>::iterator it_b[Nx][Ny];

    Grain grain;
    
    // 每个格点都是一个vector，赋值
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            // 创建结构体
            grain.key = i + 2; 
            grain.phi = 0.1;
            grid[i][j].push_back(grain);
            grid_b[i][j].push_back(grain);

            grain.key = i;
            grain.phi = 0.2;
            grid[i][j].push_back(grain);
            grid_b[i][j].push_back(grain);

            grain.key = i + 1;
            grain.phi = 0.7;
            grid[i][j].push_back(grain);
            grid_b[i][j].push_back(grain);

            it[i][j] = grid[i][j].begin();
            it_b[i][j] = grid_b[i][j].begin();
        }
    }

    // 推荐！ 通过指针访问，即使超出边界也不会崩溃，值为0。
    // cout << "Iterator access: " << endl;
    // cout << (it[185][88] + 1)->key << endl;
    // cout << (it[188][88])->phi << endl;
    // it[188][88]->phi = 6.66 + (it[185][88] + 1)->key;
    // cout << (it[188][88])->phi << endl;
    // cout << grid[188][88].at(0).phi << endl;
    // cout << endl;

    // 不推荐！通过at()访问，如果超出边界则会程序崩溃，另外，使用下标也不会崩溃，但是用指针更好
    // cout << "at access: " << endl;
    // cout << grid[185][88].at(0).key << endl;
    // cout << grid[185][88].at(0).phi << endl;
    // cout << endl;

    // 访问(185, 88)的所有结构体元素
    for (it[185][88] = grid[185][88].begin(), it_b[185][88] = grid_b[185][88].begin(); it[185][88] != grid[185][88].end(); it[185][88]++, it_b[185][88]++)
    {
        cout << it[185][88]->key << " " << it[185][88]->phi << "grid" << endl;
        cout << it_b[185][88]->key << " " << it_b[185][88]->phi << "grid_b" << endl;
    }
    cout << endl;

    // 添加结构体
    grain = {88, 0.0};
    grid[185][88].push_back(grain);
    grain = {52, 1e-7};
    grid[185][88].insert(grid[185][88].begin() + 1, grain);

    // 结构体vector排序
    sort(grid[185][88].begin(), grid[185][88].end(), LessSort);

    // 打印结果
    for (it[185][88] = grid[185][88].begin(); it[185][88] != grid[185][88].end(); it[185][88]++)
    {
        cout << it[185][88]->key << " " << it[185][88]->phi << endl;
    }
    cout << endl;

    // 除了第一个最大的数，其他的都删除。
    int tem = 0;
    for (it[185][88] = grid[185][88].begin(); it[185][88] != grid[185][88].end(); it[185][88]++)
    {
        if (tem > 0)
        {
            grid[185][88].erase(it[185][88]);
            it[185][88]--;
        }
        tem++;
    }
    cout << "delete!" << endl;

    // 打印结果
    for (it[185][88] = grid[185][88].begin(); it[185][88] != grid[185][88].end(); it[185][88]++)
    {
        cout << it[185][88]->key << " " << it[185][88]->phi << endl;
    }
    cout << endl;
    cout << "length: " << grid[185][88].size() << endl;


    double temppp = 0.17;
    int sizeee = 2;
    cout << temppp / sizeee << endl;
    cout << 0.1 << endl;

    system("pause");
}
