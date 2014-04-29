#include "graphics.h"
#include <iostream>
#include <fstream>
using namespace std;
void form3dEntity(const double* surface, int I, int J, float hx, float hy)
{
    ofstream entity;
    entity.open("entity.txt", ios::out | ios::binary);
    entity << static_cast<unsigned char>(I*J);
    entity << static_cast<unsigned char>((I-1)*(J-1));
    for (int j = 0; j<J; j++)
        for (int i = 0; i<I; i++) {
            // координаты вершины
            entity << static_cast<unsigned char>(hx*i)
                   << static_cast<unsigned char>(hy*j)
                   << static_cast<unsigned char>(surface[j*I+i]);
            // координаты нормали вершины
            entity << static_cast<unsigned char>(0)
                   << static_cast<unsigned char>(0)
                   << static_cast<unsigned char>(0);
            // текстурные координаты вершины
            entity << static_cast<unsigned char>(0)
                   << static_cast<unsigned char>(0);
        }


    for (int k = 0; k<(I-1)*J; k++) {
        // верхний треугольник
        entity << static_cast<unsigned char>(k)
               << static_cast<unsigned char>(k+I)
               << static_cast<unsigned char>(k+1)
               << static_cast<unsigned char>(0);
        // нижний треугольник
        entity << static_cast<unsigned char>(k+1)
               << static_cast<unsigned char>(k+1+I)
               << static_cast<unsigned char>(k+I)
               << static_cast<unsigned char>(0);
    }

    entity.close();
}
