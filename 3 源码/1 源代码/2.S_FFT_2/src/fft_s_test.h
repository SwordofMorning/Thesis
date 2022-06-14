#pragma once
#include <time.h>
#include "fft_s_1.h"
#include "fft_s_2.h"
#include "fft_s_4.h"

/*
*   FFT Serial 1 Test File
* 
*   文件提供两个主要的函数：
*       1. fft_s_n_test()，用于检查fft是否正确
*       2. fft_s_n_t2f()，用于检查fft的变换效率
* 
*   File provides 2 funcs:
*       1. fft_s_n_test(), to make sure the correctness of fft
*       2. fft_s_n_t2f(), statistics the runtime of fft
* 
*	1. naive fft/ifft
*	2. cooley_tukey/cooley_tukey_r
*	3. 1d/2d
*	4. r2c/c2r
*	5. convolve/correlate
*   6. radix-2/3/5
* 
*/

// 采样率
const double SampleRate = 10e3;
// 采样数量
const int SampleNums = pow(2, 10);
// radix3采样数量
const int SampleNumsRaidx3 = pow(3, 7);
// radix5采样数量
const int SampleNumsRaidx5 = pow(5, 4);
// 采样频率下界
const int SampleLowerBound = 50;
// 采样频率上界
const int SampleUpperBound = 950;

// 叠加的信号数
const int MixSigNums = 5;
// 叠加信号数组，存放信号频率
std::vector<double> MixSigVec;

// 模拟信号组数
const int SigVecs = 100;
// 模拟信号组
std::vector<std::vector<complex_t<double>>> SigVec;
// 模拟信号组 r2c使用
std::vector<std::vector<double>> SigVecReal;
// 模拟信号组 Radix3使用
std::vector<std::vector<complex_t<double>>> SigVecRadix3;
// 模拟信号组 Radix5使用
std::vector<std::vector<complex_t<double>>> SigVecRadix5;

/*  
    ==========================
    ========== Func ==========
    ==========================
*/

// fft_s_1 功能测试
void fft_s_1_test()
{
    test_convolve_fft_1d();
    test_convolve_fft_2d();
    test_convolve_fft_3d();
    size_t total_size = 8;
#if 0
    for (size_t size = 2; size <= total_size; size *= 2) {
        std::vector<complex_t<double>> t_seq;
        std::vector<complex_t<double>> f_seq;
        std::vector<complex_t<double>> t_seq_r;

        std::vector<complex_t<double>> seq_fwd;
        std::vector<complex_t<double>> seq_bwd;
        t_seq.resize(size);
        rand_vec(t_seq);
        copy_vec(t_seq, seq_fwd);

        fft_naive(t_seq, f_seq, size);
        ifft_naive(f_seq, t_seq_r, size);

        //fft_cooley_tukey(seq_fwd.data() ,size);
        fft_cooley_tukey_r(seq_fwd.data(), size);
        int err_cnt = valid_vector(f_seq, seq_fwd);

        copy_vec(f_seq, seq_bwd);
        ifft_cooley_tukey_r(seq_bwd.data(), size);
        int ierr_cnt = valid_vector(t_seq_r, seq_bwd);
        std::cout << "length:" << size << ", fwd valid:" << ((err_cnt == 0) ? "y" : "n") <<
            ", bwd valid:" << ((ierr_cnt == 0) ? "y" : "n") << std::endl;
        //dump_vector(t_seq);
        //dump_vector(t_seq_r);
        std::cout << "---------------------------------------" << std::endl;
    }
#endif
#if 1
    for (size_t size = 2; size <= total_size; size *= 2) {
        std::vector<complex_t<double>> t_seq;
        std::vector<complex_t<double>> f_seq;
        std::vector<complex_t<double>> t_seq_r;

        std::vector<double> seq_fwd_real;
        std::vector<complex_t<double>> seq_fwd;
        std::vector<complex_t<double>> seq_bwd;
        std::vector<double> seq_bwd_real;
        seq_fwd_real.resize(size);
        rand_vec(seq_fwd_real);
        for (size_t ii = 0; ii < size; ii++) {
            t_seq.push_back(complex_t<double>(seq_fwd_real[ii], (double)0));
        }

        fft_naive(t_seq, f_seq, size);
        ifft_naive(f_seq, t_seq_r, size);

        //fft_cooley_tukey(seq_fwd.data() ,size);
        seq_fwd.resize(size);
        fft_r2c(seq_fwd_real.data(), seq_fwd.data(), size);
        int err_cnt = valid_vector(f_seq, seq_fwd);

        //copy_vec(f_seq,seq_bwd);
        //ifft_cooley_tukey(seq_bwd.data(), size);
        seq_bwd_real.resize(size);
        ifft_c2r(f_seq.data(), seq_bwd_real.data(), size);
        for (size_t ii = 0; ii < size; ii++) {
            seq_bwd.push_back(complex_t<double>(seq_bwd_real[ii], (double)0));
        }
        int ierr_cnt = valid_vector(t_seq_r, seq_bwd);
        std::cout << "length:" << size << ", r2c fwd valid:" << ((err_cnt == 0) ? "y" : "n") <<
            ", c2r bwd valid:" << ((ierr_cnt == 0) ? "y" : "n") << std::endl;
        //dump_vector(t_seq);
        //dump_vector(f_seq);
        //dump_vector(seq_fwd);
        //dump_vector(t_seq_r);
        std::cout << "---------------------------------------" << std::endl;
    }
    for (size_t size = 2; size <= total_size; size *= 2) {
        size_t seq_w = size;
        size_t seq_h = size;
        std::vector<double> t_seq;
        std::vector<complex_t<double>>  f_seq;
        std::vector<double> t_seq_r;
        std::vector<complex_t<double>> t_seq_r_ex;

        std::vector<complex_t<double>>  t_seq_2;

        t_seq.resize(seq_w * seq_h);
        f_seq.resize(seq_w * seq_h);
        t_seq_r.resize(seq_w * seq_h);
        rand_vec(t_seq);
        for (size_t i = 0; i < seq_w * seq_h; i++) {
            t_seq_2.emplace_back(t_seq[i], (double)0);
        }
        fft2d_r2c(t_seq.data(), f_seq.data(), seq_w, seq_h);
        ifft2d_c2r(f_seq.data(), t_seq_r.data(), seq_w, seq_h);
        for (size_t i = 0; i < seq_w * seq_h; i++) {
            t_seq_r_ex.emplace_back(t_seq_r[i], (double)0);
        }

        fft_2d(t_seq_2.data(), seq_w, seq_h);
        int err_cnt = valid_vector(t_seq_2, f_seq);
        ifft_2d(t_seq_2.data(), seq_w, seq_h);

        int ierr_cnt = valid_vector(t_seq_2, t_seq_r_ex);

        std::cout << "length:" << seq_h << "x" << seq_w << ", r2c 2d fwd valid:" << ((err_cnt == 0) ? "y" : "n") <<
            ", c2r 2d bwd valid:" << ((ierr_cnt == 0) ? "y" : "n") << std::endl;

        std::cout << "---------------------------------------" << std::endl;
    }
#endif
}

// 随机数
int RandomNum(const int& Lowerlimit, const int& UpperLimit)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(Lowerlimit, UpperLimit);
    int num = dis(gen);
    return num;
}

// 构建模拟信号
void ConstructSignal()
{
    // 信号组
    for (int idxSigVec = 0; idxSigVec < SigVecs; ++idxSigVec)
    {
        // 清零初始化
        MixSigVec.clear();
        MixSigVec.shrink_to_fit();

        // 构造叠加信号频率
        for (int idxMixSig = 0; idxMixSig < MixSigNums; ++idxMixSig)
        {
            MixSigVec.push_back(RandomNum(SampleLowerBound, SampleUpperBound));
        }

        // 构造当前模拟信号
        std::vector<complex_t<double>> sig;
        std::vector<double> sigReal;
        
        for (int idxSig = 0; idxSig < SampleNums; ++idxSig)
        {
            double real{ 0 };

            // 叠加
            for (int idxMixSig = 0; idxMixSig < MixSigNums; ++idxMixSig)
            {
                real += sin(2 * C_PI * MixSigVec[idxMixSig] * idxSig) / SampleRate;
            }

            complex_t<double> sigVal{ real, 0 };

            sig.push_back(sigVal);
            sigReal.push_back(real);
        }

        SigVec.push_back(sig);
        SigVecReal.push_back(sigReal);
    }
}

// 构建模拟信号 Radix3
void ConstructSignalRadix3()
{
    // 信号组
    for (int idxSigVec = 0; idxSigVec < SigVecs; ++idxSigVec)
    {
        // 清零初始化
        MixSigVec.clear();
        MixSigVec.shrink_to_fit();

        // 构造叠加信号频率
        for (int idxMixSig = 0; idxMixSig < MixSigNums; ++idxMixSig)
        {
            MixSigVec.push_back(RandomNum(SampleLowerBound, SampleUpperBound));
        }

        // 构造当前模拟信号
        std::vector<complex_t<double>> sig;

        for (int idxSig = 0; idxSig < SampleNumsRaidx3; ++idxSig)
        {
            double real{ 0 };

            // 叠加
            for (int idxMixSig = 0; idxMixSig < MixSigNums; ++idxMixSig)
            {
                real += sin(2 * C_PI * MixSigVec[idxMixSig] * idxSig) / SampleRate;
            }

            complex_t<double> sigVal{ real, 0 };

            sig.push_back(sigVal);
        }

        SigVecRadix3.push_back(sig);
    }
}

// 构建模拟信号 Radix5
void ConstructSignalRadix5()
{
    // 信号组
    for (int idxSigVec = 0; idxSigVec < SigVecs; ++idxSigVec)
    {
        // 清零初始化
        MixSigVec.clear();
        MixSigVec.shrink_to_fit();

        // 构造叠加信号频率
        for (int idxMixSig = 0; idxMixSig < MixSigNums; ++idxMixSig)
        {
            MixSigVec.push_back(RandomNum(SampleLowerBound, SampleUpperBound));
        }

        // 构造当前模拟信号
        std::vector<complex_t<double>> sig;

        for (int idxSig = 0; idxSig < SampleNumsRaidx5; ++idxSig)
        {
            double real{ 0 };

            // 叠加
            for (int idxMixSig = 0; idxMixSig < MixSigNums; ++idxMixSig)
            {
                real += sin(2 * C_PI * MixSigVec[idxMixSig] * idxSig) / SampleRate;
            }

            complex_t<double> sigVal{ real, 0 };

            sig.push_back(sigVal);
        }

        SigVecRadix5.push_back(sig);
    }
}

// fft_s_1 时域2频域
void fft_s_1_t2f()
{
    clock_t start, end;
    /* ===== naive fft ===== */

    std::cout << "naive fft" << std::endl;
    start = clock();

    for (int idxSigVec = 0; idxSigVec < SigVecs; ++idxSigVec)
    {
        std::vector<complex_t<double>> output;
        fft_naive(SigVec[idxSigVec], output, SampleNums);
    }

    end = clock();
    std::cout << "time = " << double(end - start) << std::endl;
    //std::cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;
    std::cout << "==========" << std::endl;

    /* ===== cooley-tukey ===== */
    
    // 创建备份
    std::vector<std::vector<complex_t<double>>> ctCopy{ SigVec };

    std::cout << "cooley-tukey" << std::endl;
    start = clock();

    for (int idxSigVec = 0; idxSigVec < SigVecs; ++idxSigVec)
    {
        std::vector<complex_t<double>> output;
        complex_t<double>* ptr = &ctCopy[idxSigVec][0];
        fft_cooley_tukey(ptr, SampleNums);
    }

    end = clock();
    std::cout << "time = " << double(end - start) << std::endl;
    std::cout << "==========" << std::endl;

    /* ===== r2c ===== */
    
    std::cout << "r2c" << std::endl;
    start = clock();

    for (int idxSigVec = 0; idxSigVec < SigVecs; ++idxSigVec)
    {
        complex_t<double>* output = new complex_t<double>[SampleNums];
        auto ptr = &SigVecReal[idxSigVec][0];
        fft_r2c(ptr, output, SampleNums);
    }

    end = clock();
    std::cout << "time = " << double(end - start) << std::endl;
    std::cout << "==========" << std::endl;

    return;
}

// fft_s_2 功能测试
void fft_s_2_test()
{

}

// fft_s_2 时域2频域
void fft_s_2_t2f()
{
    clock_t start, end;

    /* ===== radix-2 ===== */
    // 创建备份
    std::vector<std::vector<complex_t<double>>> rdx2{ SigVec };

    std::cout << "radix-2" << std::endl;
    start = clock();

    for (int idxSigVec = 0; idxSigVec < SigVecs; ++idxSigVec)
    {
        radix_2(rdx2[idxSigVec]);
    }

    end = clock();
    std::cout << "time = " << double(end - start) << std::endl;
    std::cout << "==========" << std::endl;
    rdx2.clear();
    rdx2.shrink_to_fit();

    /* ===== radix-3 ===== */
    // 创建备份
    std::vector<std::vector<complex_t<double>>> rdx3{ SigVecRadix3 };

    std::cout << "radix-3" << std::endl;
    start = clock();

    for (int idxSigVec = 0; idxSigVec < SigVecs; ++idxSigVec)
    {
        radix_3(rdx3[idxSigVec]);
    }

    end = clock();
    std::cout << "time = " << double(end - start) << std::endl;
    std::cout << "==========" << std::endl;
    rdx3.clear();
    rdx3.shrink_to_fit();

    /* ===== radix-3 ===== */
    // 创建备份
    std::vector<std::vector<complex_t<double>>> rdx5{ SigVecRadix5 };

    std::cout << "radix-5" << std::endl;
    start = clock();

    for (int idxSigVec = 0; idxSigVec < SigVecs; ++idxSigVec)
    {
        radix_5(rdx5[idxSigVec]);
    }

    end = clock();
    std::cout << "time = " << double(end - start) << std::endl;
    std::cout << "==========" << std::endl;
    rdx5.clear();
    rdx5.shrink_to_fit();
}

// fft_s_4 时域2频域
void fft_s_4_t2f()
{
    clock_t start, end;

    /* ===== srfft-naive ===== */
    std::cout << "srfft-naive" << std::endl;
    start = clock();

    for (int idxSigVec = 0; idxSigVec < SigVecs; ++idxSigVec)
    {
        srfft1d_naive(SigVec[idxSigVec]);
    }

    end = clock();
    std::cout << "time = " << double(end - start) << std::endl;
    std::cout << "==========" << std::endl;

    /* ===== srfft-radix2 ===== */
    std::cout << "srfft-radix2" << std::endl;
    start = clock();

    for (int idxSigVec = 0; idxSigVec < SigVecs; ++idxSigVec)
    {
        srfft1d_radix2(SigVec[idxSigVec]);
    }

    end = clock();
    std::cout << "time = " << double(end - start) << std::endl;
    std::cout << "==========" << std::endl;
    
}