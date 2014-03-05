#ifndef B_AREA
#define B_AREA

class BArea{
public:
    /**
     * @brief BArea
     * @param sizeX размер области по x
     * @param sizeY размер области по y
     * @param destTime
     * @param I количество узлов по x
     * @param J количество узлов по y
     * @param T количество разбиений по времени. определяет шаг по времени dt
     */
    BArea(double sizeX, double sizeY, double destTime, int I, int J, double T);

    double answer(double x, double y, double t);

    double* H;
    int I, J, T;
    double hx, hy, dt;

    double V(double x, double y, double t);
};


#endif
