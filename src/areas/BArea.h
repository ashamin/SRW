#ifndef B_AREA
#define B_AREA

class BArea{
public:

	BArea(double sizeX, double sizeY, int I, int J);

	double answer(double x, double y, double t);

	double* H;
	int I, J;
	double hx, hy;

	double V(double x, double y, double t);
};


#endif