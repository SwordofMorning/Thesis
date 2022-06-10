#pragma once

class Complex;

class InputImage {
public:

    InputImage(const char* filename);
    int get_width() const;
    int get_height() const;

    Complex* get_image_data() const;

    void save_image_data(const char* filename, Complex* d, int w, int h);
    
    void save_image_data_real(const char* filename, Complex* d, int w, int h);

private:
    int w;
    int h;
    Complex* data;
};
