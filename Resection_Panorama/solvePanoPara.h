#ifndef SOLVEPANOPARA_H
#define SOLVEPANOPARA_H

#include <vector>

struct panoPara
{
	double xs, ys, zs;
	double alpha, phi, belta;
	double mean0[6];
	panoPara& operator=(panoPara &pp)
	{
		xs = pp.xs;
		ys = pp.ys;
		zs = pp.zs;
		alpha = pp.alpha;
		phi = pp.phi;
		belta = pp.belta;

		return *this;
	}
};

struct pointData
{
	double px, py;
	double x, y, z;
};

class SPP_S
{
public:
	SPP_S();
	SPP_S(int w, int h);

public:
	int setImgSize(int w, int h);
	int solvePanoParameter(panoPara &pp, std::vector<pointData> pd);

private:
	int imgWidth;
	int imgHeight;
	double dpi;

private:
	int computeCoefficient(panoPara pp, pointData pt, double coe[][6]);
};
#endif