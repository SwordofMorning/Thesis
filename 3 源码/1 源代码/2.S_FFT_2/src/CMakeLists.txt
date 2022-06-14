SET(
    SRC_LIST
    complex.h
    fft_s_1.h
    fft_s_2.h
    fft_s_3.h
    fft_s_4.h
    fft_s_test.h
    main.cpp
)

# 更改二进制保存目录
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR/bin})

ADD_EXECUTABLE(sfft ${SRC_LIST})