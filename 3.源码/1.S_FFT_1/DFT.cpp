#include "DFT.h"

// 构造函数
DFT::DFT(int tempSize)
{
    temp = new Comp[tempSize];
}

// 析构函数
DFT::~DFT()
{
    delete []temp;
}

void DFT::operator()(Comp* f, int n, int rev)
{
    if (n == 1) return;

    // 拷贝准备拆分奇偶
    for (int i = 0; i < n; ++i) temp[i] = f[i];

    // 偶数放左边，奇数放右边
    for (int i = 0; i < n; ++i)
    {
        if (i & 1)
            f[n / 2 + i / 2] = temp[i];
        else
            f[i / 2] = temp[i];
    }

    Comp *g = f, *h = f + n / 2;

    // 递归DFT
    this->operator()(g, n / 2, rev);
    this->operator()(h, n / 2, rev);

    Comp cur{1, 0}, step{cos(2 * M_PI / n), sin(2 * M_PI * rev / n)};

    for (int k = 0; k < n / 2; ++k)
    {
        temp[k] = g[k] + cur * h[k];
        temp[k + n / 2] = g[k] - cur * h[k];
        cur *= step; 
    }

    for (int i = 0; i < n; ++i) f[i] = temp[i];
}