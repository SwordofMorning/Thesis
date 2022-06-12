#pragma once
/*
*   FFT 单线程 4 | FFT Signal Threads
*   此文件实现算法如下：
*	1. Split-Radix，要求其长度 % 4 == 0
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
    fft1d_naive
*/
template <typename T>
std::vector<complex_t<T>> fft1d_naive(std::vector<complex_t<T>> vec)
{
    int length = vec.size();

    std::vector<complex_t<T>> fOut;
    fOut.resize(length);

    for (int k = 0; k < length; ++k)
    {
        complex_t<T> fk{ 0, 0 };
        for (int t = 0; t < length; ++t)
        {
            fk += vec[t] * omega_nkx<T>(length, k, t);
        }
        fOut[k] = fk;
    }
    return fOut;
}

/*
    split-Radix
*/
template <typename T>
std::vector<complex_t<T>> srfft1d_naive(const std::vector<complex_t<T>> vec)
{
    assert(vec.size() % 4 == 0);

    int length = vec.size();

    // init vec_2 = vec[0:length:2]
    std::vector<complex_t<T>> vec_2;
    for (int i = 0; i < length; i += 2)
    {
        vec_2.push_back(vec[i]);
    }

    // init vec_4p1 = vec[1:length:4]
    std::vector<complex_t<T>> vec_4p1;
    for (int i = 1; i < length; i += 4)
    {
        vec_4p1.push_back(vec[i]);
    }

    // init vec_4n1 = vec[3:ltngth-1:4]
    std::vector<complex_t<T>> vec_4n1;
    vec_4n1.push_back(vec[length - 1]);
    for (int i = 3; i < length - 1; i += 4)
    {
        vec_4n1.push_back(vec[i]);
    }
    
    auto fvec_2 = fft1d_naive(vec_2);
    auto fvec_4p1 = fft1d_naive(vec_4p1);
    auto fvec_4n1 = fft1d_naive(vec_4n1);
    
    std::vector<complex_t<T>> fvec;
    fvec.resize(length);

    for (int k = 0; k < length / 4; ++k)
    {
        fvec[k] =               fvec_2[k] + ((omega_nk<T>(length, k) * fvec_4p1[k] + omega_nk<T>(length, -1 * k)) * fvec_4n1[k]);
        fvec[k + length / 2] =  fvec_2[k] - ((omega_nk<T>(length, k) * fvec_4p1[k] + omega_nk<T>(length, -1 * k)) * fvec_4n1[k]);

        fvec[k + (length / 4) * 1] = fvec_2[k + length / 4] - complex_t<T>(0, 2) * (omega_nk<T>(length, k) * fvec_4p1[k] - omega_nk<T>(length, -1 * k) * fvec_4p1[k]);
        fvec[k + (length / 4) * 3] = fvec_2[k + length / 4] + complex_t<T>(0, 2) * (omega_nk<T>(length, k) * fvec_4p1[k] - omega_nk<T>(length, -1 * k) * fvec_4p1[k]);
    }

    return fvec;
}

/*
    radix-2 for srfft
*/
template<typename T>
std::vector<complex_t<T>> radix_2_4srfft(std::vector<complex_t<T>> vec)
{
    assert(is_pow_of_2(vec.size()));

    for (int i = 0; i < (vec.size() / 2); ++i)
    {
        vec[i + vec.size() / 2] = vec[i] - omega_func<T>(vec.size(), i) * vec[i + vec.size() / 2];
        vec[i] = vec[i] * 2 - vec[i + vec.size() / 2];
    }

    return vec;
}

/*
    srfft-cooley-tukey
*/
template <typename T>
std::vector<complex_t<T>> srfft1d_radix2(const std::vector<complex_t<T>> vec)
{
    assert(vec.size() % 4 == 0);

    int length = vec.size();

    // init vec_2 = vec[0:length:2]
    std::vector<complex_t<T>> vec_2;
    for (int i = 0; i < length; i += 2)
    {
        vec_2.push_back(vec[i]);
    }

    // init vec_4p1 = vec[1:length:4]
    std::vector<complex_t<T>> vec_4p1;
    for (int i = 1; i < length; i += 4)
    {
        vec_4p1.push_back(vec[i]);
    }

    // init vec_4n1 = vec[3:ltngth-1:4]
    std::vector<complex_t<T>> vec_4n1;
    vec_4n1.push_back(vec[length - 1]);
    for (int i = 3; i < length - 1; i += 4)
    {
        vec_4n1.push_back(vec[i]);
    }

    auto fvec_2 = radix_2_4srfft(vec_2);
    auto fvec_4p1 = radix_2_4srfft(vec_4p1);
    auto fvec_4n1 = radix_2_4srfft(vec_4n1);

    std::vector<complex_t<T>> fvec;
    fvec.resize(length);

    for (int k = 0; k < length / 4; ++k)
    {
        fvec[k] = fvec_2[k] + ((omega_nk<T>(length, k) * fvec_4p1[k] + omega_nk<T>(length, -1 * k)) * fvec_4n1[k]);
        fvec[k + length / 2] = fvec_2[k] - ((omega_nk<T>(length, k) * fvec_4p1[k] + omega_nk<T>(length, -1 * k)) * fvec_4n1[k]);

        fvec[k + (length / 4) * 1] = fvec_2[k + length / 4] - complex_t<T>(0, 2) * (omega_nk<T>(length, k) * fvec_4p1[k] - omega_nk<T>(length, -1 * k) * fvec_4p1[k]);
        fvec[k + (length / 4) * 3] = fvec_2[k + length / 4] + complex_t<T>(0, 2) * (omega_nk<T>(length, k) * fvec_4p1[k] - omega_nk<T>(length, -1 * k) * fvec_4p1[k]);
    }

    return fvec;
}