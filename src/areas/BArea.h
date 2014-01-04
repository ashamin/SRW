#ifndef B_AREA
#define B_AREA

class BArea{
public:

    /**
     * @param sizeX размер по X
     * @param sizeY размер по Y
     * @param destTime время
     * @param I количество разбиений по X
     * @param J количество разбиений по Y
     * @param T количество разбиений по времени
     */
	BArea(double sizeX, double sizeY, double destTime, int I, int J, double T);

	double answer(double x, double y, double t);

	double* H;
	int I, J, T;
	double hx, hy, dt;

	double V(double x, double y, double t);
};


#endif