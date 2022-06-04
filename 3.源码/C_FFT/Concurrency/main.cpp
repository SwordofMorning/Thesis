#include <iostream>
#include <fstream>
#include "DFT.h"

int main()
{
    // 采样率 10khz
    const double SAMPLE_RATE = 10e3;

    // 采样数量
    const int SAMPLE_NUMS = 1e3;

    // 模拟信号
    std::complex<double> sig[SAMPLE_NUMS];

    for (int i = 0; i < SAMPLE_NUMS; ++i)
    {
        std::complex<double> sigValue{
            // 500 + 2000
            1 * sin(2 * M_PI * 500 * i / SAMPLE_RATE) +
            1 * sin(2 * M_PI * 2000 * i / SAMPLE_RATE),
            0
        };

        sig[i] = sigValue;
    }

    std::ofstream ORI;
    ORI.open("ORI.txt");

    // 保存转换前的信号
    for (auto i : sig)
    {
        ORI << i.real() << ' ';
    }

    // 傅立叶变换
    DFT{SAMPLE_NUMS}(sig, SAMPLE_NUMS, 1);

    std::ofstream FFT;
    FFT.open("FFT.txt");

    // 保存转换后的信号
    for (auto i : sig)
    {
        FFT << sqrt(i.real() * i.real() + i.imag() * i.imag()) << ' ';
    }

    std::cout << "out" << std::endl;
    return 0;
}