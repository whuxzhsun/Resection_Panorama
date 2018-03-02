#include "stdafx.h"
#include "solvePanoPara.h"
#include <fstream>
#include <cmath>
#include "../_include/Matrix.h"
#include <iostream>

#define PI 3.1415926535897938462
#define RAD2DEG 57.29577951326

using namespace std;

SPP_S::SPP_S()
{
	ofstream outErrorLog("spp_ErrorLog.txt", ios::out);
	outErrorLog.close();
}

SPP_S::SPP_S(int w, int h)
{
	imgWidth = w;
	imgHeight = h;
	dpi = PI / imgHeight;
}

int SPP_S::setImgSize(int w, int h)
{
	imgWidth = w;
	imgHeight = h;
	dpi = PI / imgHeight;

	return 0;
}

int SPP_S::solvePanoParameter(panoPara &pp, std::vector<pointData> pd)
{
	ofstream outErrorLog("spp_ErrorLog.txt", ios::app);

	int k = 0;
	double threshould = 0;	// 循环控制阈值
	ofstream outDebug("Debug.txt", ios::out);
	outDebug.setf(ios::fixed);
	outDebug.precision(4);
	do
	{
		outDebug << "******************************* 第" << k + 1 << "次迭代 *******************************" << endl;
		int n = pd.size();
		if (n < 3)
		{
			outErrorLog << "点数小于3，不能解算全景参数！" << endl;
		}
		Matrix B(n * 2, 6);
		Matrix L(n * 2, 1);
		int i = 0, j = 0;
		for (i = 0; i < n; i++)
		{
			double coe[4][6];
			computeCoefficient(pp, pd[i], coe);

			B(i * 2 + 0, 0) = coe[0][0];	B(i * 2 + 0, 1) = coe[0][1];	B(i * 2 + 0, 2) = coe[0][2];
			B(i * 2 + 0, 3) = coe[0][3];	B(i * 2 + 0, 4) = coe[0][4];	B(i * 2 + 0, 5) = coe[0][5];
			B(i * 2 + 1, 0) = coe[1][0];	B(i * 2 + 1, 1) = coe[1][1];	B(i * 2 + 1, 2) = coe[1][2];
			B(i * 2 + 1, 3) = coe[1][3];	B(i * 2 + 1, 4) = coe[1][4];	B(i * 2 + 1, 5) = coe[1][5];
			L(i * 2 + 0, 0) = coe[2][0];	L(i * 2 + 1, 0) = coe[3][0];
		}
		Matrix Bt;
		Bt = B;
		Bt.getTranspose();
		Matrix x = (Bt * B).getInverse() * Bt * L;

		pp.xs += x(0, 0);
		pp.ys += x(1, 0);
		pp.zs += x(2, 0);
		pp.alpha += x(3, 0);
		pp.phi += x(4, 0);
		pp.belta += x(5, 0);


		/**************************************************/
		/* 输出中间结果                                   */
		/**************************************************/
		for (j = 0; j < n * 2; j++)
		{
			for (int k = 0; k < 6; k++)
			{
				outDebug << B(j, k) << "\t";
			}
			outDebug << endl;
		}
		for (j = 0; j < n * 2; j++)
		{
			outDebug << L(j, 0) << endl;
		}//*/

		Matrix v = B * x - L;

		for (j = 0; j < n; j++)
		{
			outDebug << v(j * 2 + 0, 0) << "\t" << v(j * 2 + 1, 0) << endl;
		}

		for (i = 0; i < 6; i++)
		{
			outDebug << x(i, 0) << "\t";
		}
		outDebug << endl;
		outDebug << "*************************************************************************\n" << endl;

		/*
		** 协方差因数阵、单位权中误差、外方位元素中误差
		*/
		Matrix q = (Bt * B).getInverse();
		Matrix vt = v;
		vt.getTranspose();
		double meanM0 = sqrt((vt * v)(0, 0) / (2 * n));
		for (i = 0; i < 6; i++)
		{
			pp.mean0[i] = meanM0 * sqrt(q(i, i));
		}


		threshould = (fabs(x(3, 0)) + fabs(x(4, 0)) + fabs(x(5, 0))) / 3;
		k++;
	} while (k < 15 && threshould > 0.00001);
	outDebug.close();

	if (k > 100)
	{
		outErrorLog << "相机参数解算不可靠！" << endl;
		return 1;
	}
	else
	{
		outErrorLog.close();
		return 0;
	}
}

int SPP_S::computeCoefficient(panoPara pp, pointData point, double coe[][6])
{
	double xs, ys, zs;
	xs = pp.xs;	ys = pp.ys;	zs = pp.zs;
	double ro, pt, he;
	ro = pp.alpha;	pt = pp.phi;	he = pp.belta;
	double theta, psi;
	theta = point.px * dpi - PI;
	psi = point.py * dpi;
	double x, y, z;
	x = point.x;	y = point.y;	z = point.z;

	coe[0][0] = -((cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro)) / (sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs)) + cos(pt)*sin(he)*(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs))*1.0 / pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0)) / (pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0)*1.0 / pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + 1.0);
	coe[0][1] = ((cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro)) / (sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs)) - cos(he)*cos(pt)*(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs))*1.0 / pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0)) / (pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0)*1.0 / pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + 1.0);
	coe[0][2] = ((cos(pt)*sin(ro)) / (sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs)) - sin(pt)*(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs))*1.0 / pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0)) / (pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0)*1.0 / pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + 1.0);
	coe[0][3] = -((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs)) / ((pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0)*1.0 / pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + 1.0)*(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs)));
	coe[0][4] = ((sin(pt)*sin(ro)*(z - zs) + cos(he)*cos(pt)*sin(ro)*(y - ys) + cos(pt)*sin(he)*sin(ro)*(x - xs)) / (sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs)) - (-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs))*1.0 / pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0)*(-cos(pt)*(z - zs) + cos(he)*sin(pt)*(y - ys) + sin(he)*sin(pt)*(x - xs))) / (pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0)*1.0 / pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + 1.0);
	coe[0][5] = -(((cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(x - xs) + (cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(y - ys)) / (sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs)) - (cos(he)*cos(pt)*(x - xs) - cos(pt)*sin(he)*(y - ys))*(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs))*1.0 / pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0)) / (pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0)*1.0 / pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + 1.0);

	coe[1][0] = ((((cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs))*2.0 - cos(pt)*sin(he)*(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs))*2.0)*1.0 / sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*(1.0 / 2.0)) / ((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs)) + (cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0)) / ((pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0) + 1.0);
	coe[1][1] = -((1.0 / sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*((cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs))*2.0 + cos(he)*cos(pt)*(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs))*2.0)*(1.0 / 2.0)) / ((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs)) + (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0)) / ((pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0) + 1.0);
	coe[1][2] = -(((sin(pt)*(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs))*2.0 + cos(pt)*sin(ro)*(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs))*2.0)*1.0 / sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*(1.0 / 2.0)) / ((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs)) - cos(pt)*cos(ro)*sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0)) / ((pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0) + 1.0);
	coe[1][3] = (1.0 / sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs)) + sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0)*(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs))) / ((pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0) + 1.0);
	coe[1][4] = -((((-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs))*(sin(pt)*sin(ro)*(z - zs) + cos(he)*cos(pt)*sin(ro)*(y - ys) + cos(pt)*sin(he)*sin(ro)*(x - xs))*2.0 + (sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs))*(-cos(pt)*(z - zs) + cos(he)*sin(pt)*(y - ys) + sin(he)*sin(pt)*(x - xs))*2.0)*1.0 / sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*(1.0 / 2.0)) / ((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs)) - sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0)*(cos(ro)*sin(pt)*(z - zs) + cos(he)*cos(pt)*cos(ro)*(y - ys) + cos(pt)*cos(ro)*sin(he)*(x - xs))) / ((pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0) + 1.0);
	coe[1][5] = (((((cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(x - xs) + (cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(y - ys))*(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs))*2.0 + (cos(he)*cos(pt)*(x - xs) - cos(pt)*sin(he)*(y - ys))*(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs))*2.0)*1.0 / sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*(1.0 / 2.0)) / ((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs)) + ((sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(x - xs) + (cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(y - ys))*sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0)) / ((pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0))*1.0 / pow((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs), 2.0) + 1.0);

	double tanLon = -atan((-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs)) / (sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs)));
	double tanLat = +atan(sqrt(pow(sin(pt)*(z - zs) + cos(he)*cos(pt)*(y - ys) + cos(pt)*sin(he)*(x - xs), 2.0) + pow(-(cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro))*(x - xs) + (cos(ro)*sin(he) - cos(he)*sin(pt)*sin(ro))*(y - ys) + cos(pt)*sin(ro)*(z - zs), 2.0)) / ((cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt))*(x - xs) - (sin(he)*sin(ro) + cos(he)*cos(ro)*sin(pt))*(y - ys) + cos(pt)*cos(ro)*(z - zs)));

	coe[2][0] = tanLon - theta;

	if (point.px < imgWidth / 4)
	{
		coe[2][0] = tanLon - theta - PI;
	}
	else if (point.px < imgWidth / 2)
	{
		coe[2][0] = tanLon - theta;
	}
	else if (point.px < imgWidth * 3 / 4)
	{
		coe[2][0] = tanLon - theta;
	}
	else
	{
		coe[2][0] = tanLon - theta + PI;
	}

	if (tanLat < 0)
	{
		coe[3][0] = tanLat + PI - psi;
	}
	else
	{
		coe[3][0] = tanLat - psi;
	}

	coe[2][0] = -coe[2][0];
	coe[3][0] = -coe[3][0];

	return 0;
}