#ifndef B_AREA
#define B_AREA

class BArea{
public:

    /**
     * @param sizeX ������ �� X
     * @param sizeY ������ �� Y
     * @param destTime �����
     * @param I ���������� ��������� �� X
     * @param J ���������� ��������� �� Y
     * @param T ���������� ��������� �� �������
     */
	BArea(double sizeX, double sizeY, double destTime, int I, int J, double T);

	double answer(double x, double y, double t);

	double* H;
	int I, J, T;
	double hx, hy, dt;

	double V(double x, double y, double t);
};


#endif