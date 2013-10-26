#include "BArea.h"

BArea::BArea(double sizeX, double sizeY, int I, int J){
	this->I = I;
	this->J = J;
	this->hx = sizeX / (I - 1);
	this->hy = sizeY / (J - 1);
}

double BArea::answer(double x, double y, double t){
	return x*x + y*y + t*t;
}

double BArea::V(double x, double y, double t){
	return x + y + t;
}