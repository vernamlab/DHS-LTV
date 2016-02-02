all:
    g++ -std=c++11 main.cpp general.cpp homlib.cpp ntt.cpp -o main -O3 -lntl -lgmp
