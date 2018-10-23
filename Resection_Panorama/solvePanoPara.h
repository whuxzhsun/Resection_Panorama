#ifndef SOLVEPANOPARA_H
#define SOLVEPANOPARA_H

#include <vector>

struct pointData
{
	std::string imgName;
	double px, py;
	double x, y, z;
	pointData& operator=(pointData &v)
	{
		imgName = v.imgName;
		px = v.px;	py = v.py;
		x = v.x;	y = v.y;	z = v.z;

		return *this;
	}
};

struct panoPara
{
	std::string imgName;
	double xs, ys, zs;
	double alpha, phi, belta;
	double mean0[6];
	std::vector<pointData> pixels;
	panoPara& operator=(panoPara &pp)
	{
		xs = pp.xs;
		ys = pp.ys;
		zs = pp.zs;
		alpha = pp.alpha;
		phi = pp.phi;
		belta = pp.belta;

		pixels.resize(pp.pixels.size());
		for (int i = 0; i < pp.pixels.size(); i++)
		{
			pixels[i] = pp.pixels[i];
		}

		return *this;
	}
};

class SPP_S
{
public:
	SPP_S();
	SPP_S(int w, int h);

public:
	int setImgSize(int w, int h);
	int solvePanoParameter(panoPara &pp, std::vector<pointData> pd);

	// 解算位置、姿态
	int solvePanoParameter_ceres(panoPara &pp, std::vector<pointData> pd);
	// 只解算位置
	int solvePanoParameter_ceres2(panoPara &pp, std::vector<pointData> pd);
	// 只解算姿态
	int solvePanoParameter_ceres3(panoPara &pp, std::vector<pointData> pd);

private:
	int imgWidth;
	int imgHeight;

private:
	int computeCoefficient(panoPara pp, pointData point, double coe[][6]);

};
#endif