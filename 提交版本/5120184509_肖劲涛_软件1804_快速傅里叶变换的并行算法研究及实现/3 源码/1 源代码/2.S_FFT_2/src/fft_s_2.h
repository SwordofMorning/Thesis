#pragma once
/*
*   FFT 单线程 2 | FFT Signal Threads
*   此文件实现radix分裂基算法，算法表如下：
*	1. radix-2/3/5
*	http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.50.8814
*/

#include "complex.h"
#include <vector>
#include <complex>
#include <stddef.h>
#include <utility> // std::swap in c++11
#include <assert.h>
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <string.h>
#include <functional>

/*
    旋转因子
*/
template<typename T>
auto omega_func = [](size_t total_n, size_t k)
{
    // e^( -1 * 2PI*k*n/N * i), here n is iter through each 
    T r = (T)1;
    T theta = -1 * C_2PI * k / total_n;
    return std::polar2<T>(r, theta);
};

/*
    is pow of 2
*/
bool is_pow_of_2(size_t value)
{
    return (value & (value - 1)) == 0;
}

/*
    基2，要求vec.length == pow(2, n)
*/
template<typename T>
void radix_2(std::vector<complex_t<T>>& vec)
{
    assert(is_pow_of_2(vec.size()));

    for (int i = 0; i < (vec.size() / 2); ++i)
    {
        vec[i + vec.size() / 2] = vec[i] - omega_func<T>(vec.size(), i) * vec[i + vec.size() / 2];
        vec[i] = vec[i] * 2 - vec[i + vec.size() / 2];
    }
}

/*
    2^64 = 18,446,744,073,709,551,616
    3^40 = 12,157,665,459,056,928,801
*/
bool is_pow_of_3(size_t value)
{
    return value > 0 && (12157665459056928801u) % value == 0;
}

/*
    基3，要求vec.length == pow(3, n)    
*/
template<typename T>
void radix_3(std::vector<complex_t<T>>& vec)
{
    assert(is_pow_of_3(vec.size()));

    double c1 = -1 / 2;
    double c2 = (sqrt(3) / 2) * (-1);

    /*
        z1 = omega(length, j) * vec[j + length//3]
        s1 = z1 - omega(length, 2*j)*vec[j + (length*2)//3]
        s2 = 2*z1 - s1
        s3 = s2 + vec[j]
        s4 = vec[j] + c1*s2
        s5 = s4 - c2*s1
        s6 = 2*s4 - s5
        vec[j] = s3
        vec[j+length//3] = s6
        vec[j+(2*length)//3] = s5
    */

    int length = vec.size();

    for (int i = 0; i < (length / 3); ++i)
    {
        vec[i + length / 3] = omega_func<T>(length, i) * vec[i + length / 3];
        vec[i + (length * 2) / 3] = vec[i + length / 3] - omega_func<T>(length, i) * vec[i + (length * 2) / 3];
        vec[i + length / 3] = vec[i + length / 3] * 2 - vec[i + (length * 2) / 3];
        auto tmp0 = vec[i] + vec[i + length / 3] * c1;
        vec[i] = vec[i + length / 3] + vec[i];
        vec[i + (length * 2) / 3] = tmp0 - vec[i + (length * 2) / 3] * c2;
        vec[i + length / 3] = tmp0 * 2 - vec[i + (length * 2) / 3];
    }
}

/*
    2^64 = 18,446,744,073,709,551,616
    5^27 = 7,450,580,596,923,828,125
*/
bool is_pow_of_5(size_t value)
{
    return value > 0 && (7450580596923828125u) % value == 0;
}

/*
    基5，要求vec.length == pow(5, n)
*/
template<typename T>
void radix_5(std::vector<complex_t<T>>& vec)
{
    assert(is_pow_of_5(vec.size()));

    double c1 = 1 / 4;
    double c2 = sqrt(5) / 4;
    double c3 = sqrt(
        (5 - sqrt(5)) /
        (5 + sqrt(5)));
    double c4 = (1 / 2) * sqrt(5 / 2 + sqrt(5) / 2);

    int length = vec.size();

    for (int i = 0; i < length / 5; ++i)
    {
        auto z0 = vec[i];
        auto z1 = omega_func<T>(length, i) + vec[i + length / 5];
        auto z2 = omega_func<T>(length, 2 * i) * vec[i + (length * 2) / 5];
        auto s1 = z1 - omega_func<T>(length, 4 * i) * vec[i + (length * 4) / 5];
        auto s2 = z1 * 2 - s1;
        auto s3 = z2 - omega_func<T>(length, 3 * i) * vec[i + (length * 3) / 5];
        auto s4 = z2 * 2 - s3;
        auto s5 = s2 + s4;
        auto s6 = s2 - s4;
        auto s7 = z0 - s5 * c1;
        auto s8 = s7 - s6 * c2;
        auto s9 = s7 * 2 - s8;
        auto s10 = s1 + s3 * c3;
        auto s11 = s1 * c3 - s3;
        
        vec[i] = z0 + s5;
        auto t1 = s9 - s10 * c4;
        vec[i + (length * 4) / 5] = s9 * 2 - t1;
        auto t2 = s8 - s11 * c4;
        vec[i + (length * 3) / 5] = s8 * 2 - t2;
        vec[i + (length * 2) / 5] = t2;
        vec[i + (length * 1) / 5] = t1;
    }
}