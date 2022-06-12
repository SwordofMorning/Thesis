#include "fft_s_1.h"
#include "fft_s_2.h"
#include "fft_s_4.h"
#include "fft_s_test.h"

void FFT_Serial()
{
    ConstructSignal();
    fft_s_1_t2f();

    ConstructSignalRadix3();
    ConstructSignalRadix5();
    fft_s_2_t2f();

    fft_s_4_t2f();
}

int main(int argc, char* argv[]) 
{
    FFT_Serial();
    return 0;
}