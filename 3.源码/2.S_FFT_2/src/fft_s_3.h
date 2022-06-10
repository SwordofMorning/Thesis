#pragma once
/*
*   FFT 单线程 3 | FFT Signal Threads
*   此文件实现算法如下：
*	1. 互质因子算法
*	PFA https://enacademic.com/dic.nsf/enwiki/151599
*	https://pdfs.semanticscholar.org/18e9/67f8f17ef30e6ab8f77b0f6fe56b0af4abd4.pdf
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
#include <map>

template <typename T>
auto omega_nk = [](size_t total_n, size_t k)
{
    // n is total pixel, k is current requested frequency at k.
    T theta = -1 * C_2PI * k / total_n;
    return complex_t<T>(cos(theta), sin(theta));
};

template <typename T>
auto omega_nkx = [](size_t total_n, size_t k, size_t x)
{
    // omega at input x
    T theta = -1 * C_2PI * k * x / total_n;
    return complex_t<T>(cos(theta), sin(theta));
};

/*
    fft_1d_naive, int fft_s_1.h
    直接调用fft_s_1中的fft_naive()函数即可
*/

/*
    greatest common divisor
    计算最大公因数
*/
int gcd(int p, int q)
{
    while (q)
    {
        /*
            p = q;
            q = p % q;
        */
        int x = p;
        p = q;
        q = x % q;
    }
    return p;
}

/*
    linear congruence
    线性同余方程求解
    http://gauss.math.luc.edu/greicius/Math201/Fall2012/Lectures/linear-congruences.article.pdf
    a*x ~= b (mod m)
    https://zh.wikipedia.org/wiki/%E7%BA%BF%E6%80%A7%E5%90%8C%E4%BD%99%E6%96%B9%E7%A8%8B
*/
std::vector<int> linear_congruence(int a, int b, int m)
{
    // 返回值
    std::vector<int> re;

    // 边界值检查
    if (b == 0) return re;
    assert(a > 0 && b > 0 && m > 0);

    // 计算公因数d
    int d = gcd(a, m);

    // 无解
    if (b % d != 0)
    {
        std::cout << "无解" << std::endl;
        return re;
    }

    int x = 0;
    for (; x < m; ++x)
    {
        if ((a * x - b) % m == 0)
            break;
    }

    // 一个解
    if (d == 1)
    {
        re.push_back(x);
        return re;
    }

    // 多个解
    for (int i = 0; i < d; ++i)
    {
        re.push_back(x + (m / d) * i);
    }
    return re;
}

/*
    丢番图方程
    http://volta.sdsu.edu/~amir/Bhagat18High.pdf
    need solve equation to find m1, m2
    n1*m1 + n2*m2 ~= 1 mod (n1*n2)
    -> n1*m1 ~= 1 mod n2
    -> n2*m2 ~= 1 mod n1
*/
std::pair<int, int> solve_mapping_congruence_diophantine(int map_n1, int map_n2)
{
    auto map_m1 = linear_congruence(map_n1, 1, map_n2);
    auto map_m2 = linear_congruence(map_n2, 1, map_n1);
    
    assert(!map_m1.empty() && !map_m2.empty());

    auto m1 = map_m1[0];
    auto m2 = map_m2[0];

    assert((map_n1 * m1 * map_n2 * m2) % (map_n1 * map_n2) == 1);
    
    return std::make_pair(m1, m2);
}

/*  FWD
    ruritanian correspondence mapping
    n1, n2 ro n
*/
std::vector<std::vector<int>> rcm_fwd_mapping(int map_n1, int map_n2)
{
    // need to be co-prime
    assert(gcd(map_n1, map_n2) == 1);

    // 创建一个n1 * n2大小的二维数组
    std::vector<std::vector<int>> mapping;
    mapping.resize(map_n1);
    for (int i = 0; i < map_n1; ++i)
    {
        mapping[i].resize(map_n2);
    }

    // 初始化
    for (int idx1 = 0; idx1 < map_n1; ++idx1)
    {
        for (int idx2 = 0; idx2 < map_n2; ++idx2)
        {
            mapping[idx1][idx2] = (idx1 * map_n2 + idx2 * map_n1) % (map_n1 * map_n2);
        }
    }

    return mapping;
}

/*
    BWD
    ruritanian correspondence mapping
    n to n1, n2
*/
std::vector<std::vector<int>> rcm_bwd_mapping(int map_n1, int map_n2)
{
    // need to be co-prime
    assert(gcd(map_n1, map_n2) == 1);

    auto map_m_pair = solve_mapping_congruence_diophantine(map_n1, map_n2);
    int map_m1 = map_m_pair.first;
    int map_m2 = map_m_pair.second;

    // 创建一个n1 * n2大小的二维数组
    std::vector<std::vector<int>> mapping;
    mapping.resize(map_n1);
    for (int i = 0; i < map_n1; ++i)
    {
        mapping[i].resize(map_n2);
    }

    // 初始化
    for (int n = 0; n < map_n1 * map_n2; ++n) 
    {
        int in1 = (n * map_m2) % map_n1;
        int in2 = (n * map_m1) % map_n2;
        mapping[n][0] = in1;
        mapping[n][1] = in2;
    }

    return mapping;
}

/*
    余数定理
    Fwd
    n to n1, n2
*/
std::vector<std::vector<int>> crt_fwd_mapping(int map_n1, int map_n2)
{
    // need to be co-prime
    assert(gcd(map_n1, map_n2) == 1);

    // 创建一个n1 * n2大小的二维数组
    std::vector<std::vector<int>> mapping;
    mapping.resize(map_n1);
    for (int i = 0; i < map_n1; ++i)
    {
        mapping[i].resize(map_n2);
    }

    auto map_m_pair = solve_mapping_congruence_diophantine(map_n1, map_n2);
    int map_m1 = map_m_pair.first;
    int map_m2 = map_m_pair.second;

    for (int idx1 = 0; idx1 < map_n1; idx1++)
    {
        for (int idx2 = 0; idx2 < map_n2; ++idx2)
        {
            mapping[idx1][idx2] = (idx1 * map_n2 * map_m2 + idx2 * map_n1* map_m1) % (map_n1 * map_n2);
        }
    }
    return mapping;
}

/*
    余数定理
    bwd
    n1, n2 to n
*/
std::vector<std::vector<int>> crt_bwd_mapping(int map_n1, int map_n2)
{
    // need to be co-prime
    assert(gcd(map_n1, map_n2) == 1);

    // 创建一个n1 * n2大小的二维数组
    std::vector<std::vector<int>> mapping;
    mapping.resize(map_n1);
    for (int i = 0; i < map_n1; ++i)
    {
        mapping[i].resize(map_n2);
    }

    for (int n = 0; n < map_n1 * map_n2; ++n)
    {
        int in1 = n % map_n1;
        int in2 = n % map_n2;
        mapping[n][0] = in1;
        mapping[n][1] = in2;
    }

    return mapping;
}

/*
    map_1d_to_2d
*/
template <typename T>
std::vector<std::vector<complex_t<T>>> map_1d_to_2d(std::vector<complex_t<T>> vec_1d, std::vector<std::vector<int>> mapping)
{
    // 宽高 Map's Height and Width
    int map_n1 = mapping.size();
    int map_n2 = mapping[0].size();

    assert(vec_1d.size() == map_n1 * map_n2);

    std::vector<std::vector<complex_t<T>>> vec_2d;
    vec_2d.resize(map_n1);
    for (int idx = 0; idx < map_n1; ++idx)
    {
        vec_2d[idx].resize(map_n2);
    }


    for (int idx1 = 0; idx1 < map_n1; ++idx1)
    {
        for (int idx2 = 0; idx2 < map_n2; ++idx2)
        {
            vec_2d[idx1][idx2] = vec_1d(mapping[idx1][idx2]);
        }
    }

    return vec_2d;
}



/*
    PFA
*/
template <typename T>
void pfa_fft_1d_naive_n1_n2(std::vector<complex_t<T>> vec, int n1, int n2)
{

}