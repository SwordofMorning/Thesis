#include <complex>
#include <cmath>

// 离散傅立叶变换
class DFT
{
private:
    using Comp = std::complex<double>;

    // 临时拷贝用的成员变量
    Comp* temp;

    // 虚数i
    const Comp I{0, 1};

public:
    DFT(int tempSize);

    ~DFT();
    
    void operator()(Comp* f, int n, int rev);
};