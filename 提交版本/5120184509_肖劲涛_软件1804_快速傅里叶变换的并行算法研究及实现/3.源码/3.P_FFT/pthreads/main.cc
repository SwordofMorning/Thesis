#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <thread>
#include <Windows.h>
#include <process.h>
#include "complex.h"

const float PI = 3.14159265358979f;
static const int N_THREADS = 8;

std::vector<Complex> compute_row(const std::vector<Complex> &inputArray, int dim)
{
    std::vector<Complex> dft(dim);
    for (int i = 0; i < dim; i++) 
    {
        Complex sum;
        for (int j = 0; j < dim; j++)
        {
            Complex W(cos(2 * PI * i * j / dim), -1 * sin(2 * PI * i * j / dim));
            Complex in(inputArray[j]);
            sum = sum + (W * in);
        }
        dft[i] = sum;
    }
    return dft;
}

std::vector<Complex> rcompute_row(const std::vector<Complex> &inputArray, int dim)
{
    std::vector<Complex> dft(dim);
    for (int i = 0; i < dim; i++) 
    {
        Complex sum;
        for (int j = 0; j < dim; j++)
        {
            Complex W(cos(2 * PI * i * j / dim), sin(2 * PI * i * j / dim));
            Complex in(inputArray[j]);
            sum = sum + (W * in);
        }
        dft[i] = sum * (1 / float(dim));
    }
    return dft;
}

struct Arg
{
    const std::vector<std::vector<Complex>>& inputArray;
    std::vector<std::vector<Complex>>& newArray;
    int start;
    int stop;
    int dim;
};

unsigned _stdcall threaded_dft(void* arglist)
{
    Arg* p = (Arg*)arglist;
    for (int i = p->start; i < p->stop; i++)
    {
        p->newArray[i] = compute_row(p->inputArray[i], p->dim);
    }
    return 1;
}

void threaded_rdft(const std::vector<std::vector<Complex>> &inputArray, std::vector<std::vector<Complex>> &newArray, int start, int stop, int dim)
{
    for (int i = start; i < stop; i++)
    {
        newArray[i] = rcompute_row(inputArray[i], dim);
    }
}

std::vector<std::string> read_input(std::string configFile)
{
    std::ifstream config;
    config.open(configFile);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(config, line))
    {
        if (line.front() != '#' && !line.empty())
        {
            lines.push_back(line);
        }
    }
    return lines;
}

void print_array(std::string output_file, int dim, const std::vector<std::vector<Complex>> &outArray)
{
    std::ofstream myfile;
    myfile.open(output_file);
    myfile << dim << " " << dim << std::endl;
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            myfile << outArray[i][j] << " ";
        }
        myfile << std::endl;
    }
}

std::vector<Complex> parse_string(std::string str, std::string delim)
{
    size_t pos = 0;
    std::vector<Complex> complexs;
    while ((pos = str.find(delim)) != std::string::npos)
    {
        complexs.push_back(atoi(str.substr(0, pos).c_str()));
        str.erase(0, pos + delim.length());
    }
    complexs.push_back(atoi(str.c_str()));
    return complexs;
}

void dft_row(int dim, std::vector<std::vector<Complex>> &inputArray, std::vector<std::vector<Complex>> &rowArray)
{
    int rows_per = dim / N_THREADS;
    int remainder = dim % N_THREADS;

    
    HANDLE t[N_THREADS];
    for (int i = 0; i < N_THREADS; i++)
    {
        
        t[i] = (HANDLE)_beginthreadex(nullptr, 0, threaded_dft, nullptr, 0, nullptr);
        SetThreadAffinityMask(t[i], i);
    }

    for (int i = dim - remainder; i < dim; i++)
    {
        rowArray[i] = compute_row(inputArray[i], dim);
    }

    for (auto& th : t)
        WaitForSingleObject(th, INFINITE);

    for (auto& th : t)
        CloseHandle(th);
}

void rdft_row(int dim, std::vector<std::vector<Complex>> &inputArray, std::vector<std::vector<Complex>> &rowArray)
{
    int rows_per = dim / N_THREADS;
    int remainder = dim % N_THREADS;

    std::thread t[N_THREADS];
    for (int i = 0; i < N_THREADS; i++)
    {
        t[i] = std::thread(threaded_rdft, std::ref(inputArray), std::ref(rowArray), rows_per * i, (rows_per * (i + 1)), dim);
    }

    for (int i = dim - remainder; i < dim; i++)
    {
        rowArray[i] = rcompute_row(inputArray[i], dim);
    }

    for (auto &th : t)
        th.join();
}

std::vector<std::vector<Complex>> transpose(int dim, const std::vector<std::vector<Complex>> &inputArray)
{
    std::vector<std::vector<Complex>> transposedArray(dim, std::vector<Complex>(dim));
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            transposedArray[j][i] = inputArray[i][j];
    return transposedArray;
}

void run_forward(std::vector<std::string> input_data, std::string output_file)
{
    
    int dim = parse_string(input_data[0], " ")[0].real;

    
    std::vector<std::vector<Complex>> inputArray(dim, std::vector<Complex>(dim));
    std::vector<std::vector<Complex>> rowArray(dim, std::vector<Complex>(dim));
    std::vector<std::vector<Complex>> colArray(dim, std::vector<Complex>(dim));
    for (int i = 1; i < dim + 1; i++)
    {
        inputArray[i - 1] = parse_string(input_data[i], " ");
    }

    
    dft_row(dim, inputArray, rowArray);
    colArray = transpose(dim, rowArray);
    dft_row(dim, colArray, inputArray);
    rowArray = transpose(dim, inputArray);

    
    print_array(output_file, dim, rowArray);
}

void run_reverse(std::vector<std::string> input_data, std::string output_file)
{
    
    int dim = parse_string(input_data[0], " ")[0].real;

    
    std::vector<std::vector<Complex>> inputArray(dim, std::vector<Complex>(dim));
    std::vector<std::vector<Complex>> rowArray(dim, std::vector<Complex>(dim));
    std::vector<std::vector<Complex>> colArray(dim, std::vector<Complex>(dim));
    for (int i = 1; i < dim + 1; i++)
    {
        inputArray[i - 1] = parse_string(input_data[i], " ");
    }

    
    rdft_row(dim, inputArray, rowArray);
    colArray = transpose(dim, rowArray);
    rdft_row(dim, colArray, inputArray);
    rowArray = transpose(dim, inputArray);

    
    print_array(output_file, dim, rowArray);
}

int main(int argc, char **argv)
{
    
    
    
    std::vector<std::string> input_data = read_input(argv[2]);
    std::string dir(argv[1]);
    std::string output_file(argv[3]);

    

    std::clock_t start, end;

    start = clock();

    if (dir == "forward")
    {
        run_forward(input_data, output_file);

        std::cout << "Forward: " << N_THREADS << std::endl;
    }
    else if (dir == "reverse")
    {
        run_reverse(input_data, output_file);
    }

    end = clock();

    std::cout << double(end - start) << std::endl;

    return 1;
}
