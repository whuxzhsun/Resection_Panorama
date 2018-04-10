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

		Matrix P(2 * n, n * 2);
		for (i = 0; i < n * 2; i++)
		{
			for (j = 0; j < n * 2; j++)
			{
				if (i == j)
				{
					P(i, j) = sqrt(pow(pd[i / 2].x - pp.xs, 2) + pow(pd[i / 2].y - pp.ys, 2) + pow(pd[i /2].z - pp.zs, 2));
				}
				else
				{
					P(i, j) = 0;
				}
			}
		}

		Matrix Bt;
		Bt = B;
		Bt.getTranspose();
		Matrix x = (Bt * P * B).getInverse() * Bt * P * L;

		pp.xs += x(0, 0);
		pp.ys += x(1, 0);
		pp.zs += x(2, 0);
		pp.alpha	+= x(3, 0);
		pp.phi		+= x(4, 0);
		pp.belta	+= x(5, 0);

		// 输出中间过程以供查看
		outDebug.setf(ios::fixed);
		for (j = 0; j < n * 2; j++)
		{
			for (int k = 0; k < 6; k++)
			{
				outDebug.width(6);
				outDebug.precision(3);
				outDebug << B(j, k) << "\t";
			}
			outDebug << endl;
		}
		for (j = 0; j < n * 2; j++)
		{
			outDebug.width(6);
			outDebug.precision(3);
			outDebug << L(j, 0) << endl;
		}

		Matrix v = B * x - L;

		for (j = 0; j < n; j++)
		{
			pd[j].px += v(j * 2 + 0, 0) / dpi;
			pd[j].py += v(j * 2 + 1, 0) / dpi;
// 			cout.setf(ios::fixed);
// 			cout.width(6);
// 			cout.precision(3);
// 			cout << v(j * 2 + 0, 0) / dpi << "\t";
// 			cout.width(6);
// 			cout.precision(3);
// 			cout << v(j * 2 + 1, 0) / dpi << endl;
			
			outDebug.width(6);
			outDebug.precision(3);
			outDebug << v(j * 2 + 0, 0) / dpi << "\t";
			outDebug.width(6);
			outDebug.precision(3);
			outDebug << v(j * 2 + 1, 0) / dpi << endl;
		}
/*		cout << endl;*/

		for (i = 0; i < 6; i++)
		{
			outDebug.width(6);
			outDebug.precision(3);
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
	} while (++k < 15 && threshould > 0.00001);
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
	double ro, pt, he;
	ro = pp.alpha;	
	pt = pp.phi;	
	he = pp.belta;

	double theta, psi;
	theta = point.px * dpi - PI;
	psi = point.py * dpi;

	double r1, r2, r3, r4, r5, r6, r7, r8, r9;
	r1 = cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro);
	r2 = cos(he)*sin(pt)*sin(ro) - cos(ro)*sin(he);
	r3 =-cos(pt)*sin(ro);
	r4 = cos(pt)*sin(he);
	r5 = cos(he)*cos(pt);
	r6 = sin(pt);
	r7 = cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt);
	r8 =-sin(he)*sin(ro) - cos(he)*cos(ro)*sin(pt); 
	r9 = cos(pt)*cos(ro);

	double bx, by, bz;
	bx = point.x - pp.xs;
	by = point.y - pp.ys;
	bz = point.z - pp.zs;

	double A, B, C;
	A = r1*bx + r2*by + r3*bz;
	B = r4*bx + r5*by + r6*bz;
	C = r7*bx + r8*by + r9*bz;

	double A2, B2, C2, A2B2, SQAB;
	A2 = A * A;
	B2 = B * B;
	C2 = C * C;
	A2B2 = A2 + B2;
	SQAB = (A2B2 + C2) * sqrt(A2B2);

	double tp1, tp2, tp3, tp4, tp5, tp6, tp7, tp8, tp9, tp10, tp11, tp12, tp13;
	tp1 = sin(pt)*sin(ro);
	tp2 = cos(he)*cos(pt)*sin(ro);
	tp3 = cos(pt)*sin(he)*sin(ro);
	tp4 = cos(he)*sin(pt);
	tp5 = cos(pt);
	tp6 = sin(he)*sin(pt);
	tp7 = cos(ro)*sin(pt);
	tp8 = cos(he)*cos(pt)*cos(ro);
	tp9 = cos(pt)*cos(ro)*sin(he);
	tp10 = C * (A * (tp1*bz + tp2*by + tp3*bx) - B * (tp4*by - tp5*bz + tp6*bx));
	tp11 = A2B2 * (tp7*bz + tp8*by + tp9*bx);
	tp12 = C * (A * (r2*bx - r1*by) + B * (r5*bx - r4*by));
	tp13 = A2B2 * (-r8*bx + r7*by);

	coe[0][0] = (A * r4 - B * r1) / A2B2;
	coe[0][1] = (A * r5 - B * r2) / A2B2;
	coe[0][2] = (A * r6 - B * r3) / A2B2;
	coe[0][3] = -B * C / A2B2;
	coe[0][4] = (A * (tp4 * by - tp5 * bz + tp6 * bx) + B * (tp1 * bz + tp2 * by + tp3 * bx)) / A2B2;
	coe[0][5] = (A * (r4 * by - r5 * bx) + B * (r2 * bx - r1 * by)) / A2B2;

	coe[1][0] = (A2B2 * r7 - C * (A * r1 + B * r4)) / SQAB;
	coe[1][1] = (A2B2 * r8 - C * (A * r2 + B * r5)) / SQAB;
	coe[1][2] = (A2B2 * r9 - C * (A * r3 + B * r6)) / SQAB;
	coe[1][3] = -A / sqrt(A2B2);
	coe[1][4] = (tp10 + tp11) / SQAB;
	coe[1][5] = (tp12 + tp13) / SQAB;

	double tanLon = 0, tanLat = 0;
	tanLon = atan(A / B);
	tanLat = atan(sqrt(A2B2) / C);

	if (point.px < imgWidth / 4)
	{
		coe[2][0] = tanLon - PI - theta;
	}
	else if (point.px > imgWidth * 3 / 4)
	{
		coe[2][0] = tanLon + PI - theta;
	}
	else
	{
		coe[2][0] = tanLon - theta;
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