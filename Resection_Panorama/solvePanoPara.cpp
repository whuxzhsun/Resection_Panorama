#include "stdafx.h"
#include "solvePanoPara.h"
#include <fstream>
#include <cmath>
#include "../_include/Matrix.h"
#include <iostream>

#include "ceres/ceres.h"
#include "glog/logging.h"

#define PI 3.1415926535897938462
#define RAD2DEG 57.29577951326
#define DPI  (PI / 4096)

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
}

int SPP_S::setImgSize(int w, int h)
{
	imgWidth = w;
	imgHeight = h;

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
			pd[j].px += v(j * 2 + 0, 0) / DPI;
			pd[j].py += v(j * 2 + 1, 0) / DPI;
// 			cout.setf(ios::fixed);
// 			cout.width(6);
// 			cout.precision(3);
// 			cout << v(j * 2 + 0, 0) / DPI << "\t";
// 			cout.width(6);
// 			cout.precision(3);
// 			cout << v(j * 2 + 1, 0) / DPI << endl;
			
			outDebug.width(6);
			outDebug.precision(3);
			outDebug << v(j * 2 + 0, 0) / DPI << "\t";
			outDebug.width(6);
			outDebug.precision(3);
			outDebug << v(j * 2 + 1, 0) / DPI << endl;
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
	theta = point.px * DPI - PI;
	psi = point.py * DPI;

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
	tanLon = atan2(A, B);
	tanLat = atan(sqrt(A2B2) / C);

	coe[2][0] = tanLon - theta;
// 	if (point.px < imgWidth / 4)
// 	{
// 		coe[2][0] = tanLon - PI - theta;
// 	}
// 	else if (point.px > imgWidth * 3 / 4)
// 	{
// 		coe[2][0] = tanLon + PI - theta;
// 	}
// 	else
// 	{
// 		coe[2][0] = tanLon - theta;
// 	}

	if (coe[2][0] < -PI / 2)
	{
		coe[2][0] += PI;
	}
	else if (coe[2][0] > PI / 2)
	{
		coe[2][0] -= PI;
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

struct ResectionResidual
{
	ResectionResidual(double px, double py, double x, double y, double z)
	:_px(px), _py(py), _x(x), _y(y), _z(z)
	{}

	template <typename T>
	bool operator()(const T * const parameter, T *residual) const
	{
		T ro, pt, he;
		ro = parameter[3];
		pt = parameter[4];
		he = parameter[5];

		T r1, r2, r3, r4, r5, r6, r7, r8, r9;
		r1 = cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro);
		r2 = cos(he)*sin(pt)*sin(ro) - cos(ro)*sin(he);
		r3 = -cos(pt)*sin(ro);
		r4 = cos(pt)*sin(he);
		r5 = cos(he)*cos(pt);
		r6 = sin(pt);
		r7 = cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt);
		r8 = -sin(he)*sin(ro) - cos(he)*cos(ro)*sin(pt);
		r9 = cos(pt)*cos(ro);

		T bx, by, bz;
		bx = _x - parameter[0];
		by = _y - parameter[1];
		bz = _z - parameter[2];

		T A, B, C;
		A = r1*bx + r2*by + r3*bz;
		B = r4*bx + r5*by + r6*bz;
		C = r7*bx + r8*by + r9*bz;

		T tanLon = atan2(A, B);
		T tanLat = atan2(sqrt(A * A + B * B), C);

 		residual[0] = tanLon - T(_px);
		residual[1] = tanLat - T(_py);

// 		if (residual[0] < T(-PI / 2))
// 		{
// 			residual[0] = residual[0] + T(PI);
// 		}
// 		else if (residual[0] > T(PI / 2))
// 		{
// 			residual[0] = residual[0] - T(PI);
// 		}
// 
// 		if (tanLat < T(0))
// 		{
// 			residual[1] = tanLat + T(PI);
// 		}

		return true;
	}

private:
	const double _px, _py, _x, _y, _z;
};

int SPP_S::solvePanoParameter_ceres(panoPara &pp, std::vector<pointData> pd)
{
	double para[6] = { pp.xs, pp.ys, pp.zs, pp.alpha, pp.phi, pp.belta };

	ceres::Problem problem;
	for (int i = 0; i < pd.size(); i++)
	{
		pd[i].px = pd[i].px * DPI - PI;
		pd[i].py = pd[i].py * DPI;
		problem.AddResidualBlock(new ceres::AutoDiffCostFunction<ResectionResidual, 2, 6>
			(new ResectionResidual(pd[i].px, pd[i].py, pd[i].x, pd[i].y, pd[i].z)), 
			/*new ceres::LossFunctionWrapper(new ceres::SoftLOneLoss(1.0), ceres::TAKE_OWNERSHIP)*/
			NULL, para);
	}

	ceres::Solver::Options options;
	options.max_lbfgs_rank = 20;
	options.use_approximate_eigenvalue_bfgs_scaling = false;
	options.min_line_search_step_size = 1e-9;
	options.line_search_sufficient_function_decrease = 1e-4;
	options.max_line_search_step_contraction = 1e-3;
	options.min_line_search_step_contraction = 0.6;
	options.max_num_line_search_step_size_iterations = 20;
	options.max_num_line_search_direction_restarts = 5;
	options.line_search_sufficient_curvature_decrease = 0.9;
	options.max_line_search_step_expansion = 10.0;
	options.use_nonmonotonic_steps = false;
	options.max_consecutive_nonmonotonic_steps = 5;
	options.max_num_iterations = 50;
	options.initial_trust_region_radius = 1e4;
	options.max_trust_region_radius = 1e16;
	options.min_trust_region_radius = 1e-32;
	options.min_relative_decrease = 1e-3;	// default 1e-3
	options.min_lm_diagonal = 1e-6;
	options.max_lm_diagonal = 1e32;
	options.max_num_consecutive_invalid_steps = 5;
	options.function_tolerance = 1e-8;	// default 1e-6;
	options.gradient_tolerance = 1e-12;	// default 1e-10;
	options.parameter_tolerance = 1e-10;// defautl le-8
	options.gradient_check_relative_precision = 1e-8;	// default le-8
	options.numeric_derivative_relative_step_size = 1e-6;	// default le-6
	options.minimizer_progress_to_stdout = false;

	ceres::Solver::Summary summary;

	ceres::Solve(options, &problem, &summary);

	pp.xs = para[0];
	pp.ys = para[1];
	pp.zs = para[2];
	pp.alpha = para[3];
	pp.phi = para[4];
	pp.belta = para[5];

	return 0;
}

struct ResectionResidual2
{
	ResectionResidual2(double px, double py, double x, double y, double z, double r, double p, double h)
	:_px(px), _py(py), _x(x), _y(y), _z(z), _r(r), _p(p), _h(h)
	{}

	template <typename T>
	bool operator()(const T * const parameter_0, /*const T * const parameter_1,*/ T *residual) const
	{
		T ro, pt, he;
		ro = T(_r);
		pt = T(_p);
		he = T(_h);

		T r1, r2, r3, r4, r5, r6, r7, r8, r9;
		r1 = cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro);
		r2 = cos(he)*sin(pt)*sin(ro) - cos(ro)*sin(he);
		r3 = -cos(pt)*sin(ro);
		r4 = cos(pt)*sin(he);
		r5 = cos(he)*cos(pt);
		r6 = sin(pt);
		r7 = cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt);
		r8 = -sin(he)*sin(ro) - cos(he)*cos(ro)*sin(pt);
		r9 = cos(pt)*cos(ro);

		T bx, by, bz;
		bx = _x - parameter_0[0];
		by = _y - parameter_0[1];
		bz = _z - parameter_0[2];

		T A, B, C;
		A = r1*bx + r2*by + r3*bz;
		B = r4*bx + r5*by + r6*bz;
		C = r7*bx + r8*by + r9*bz;

		T tanLon = atan2(A, B);
		T tanLat = atan2(sqrt(A * A + B * B), C);

		residual[0] = tanLon - T(_px);
		residual[1] = tanLat - T(_py);

		return true;
	}

private:
	const double _px, _py, _x, _y, _z;
	const double _r, _p, _h;
};

// 只解算xyz
int SPP_S::solvePanoParameter_ceres2(panoPara &pp, std::vector<pointData> pd)
{
	double para_0[3] = {pp.xs, pp.ys, pp.zs};
	double para_1[3] = {pp.alpha, pp.phi, pp.belta};

	ceres::Problem problem;
	for (int i = 0; i < pd.size(); i++)
	{
		pd[i].px = pd[i].px * DPI - PI;
		pd[i].py = pd[i].py * DPI;

		problem.AddResidualBlock(new ceres::AutoDiffCostFunction<ResectionResidual2, 2, 3>(new ResectionResidual2(pd[i].px, pd[i].py, pd[i].x, pd[i].y, pd[i].z, pp.alpha, pp.phi, pp.belta)),
			NULL, para_0/*, para_1*/);
	}

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	ceres::Solve(options, &problem, &summary);

	pp.xs = para_0[0];
	pp.ys = para_0[1];
	pp.zs = para_0[2];
	pp.alpha = para_1[0];
	pp.phi = para_1[1];
	pp.belta = para_1[2];

	return 0;
}

struct ResectionResidual3
{
	ResectionResidual3(double px, double py, double x, double y, double z, double r, double p, double h)
	:_px(px), _py(py), _x(x), _y(y), _z(z), _xs(r), _ys(p), _zs(h)
	{}

	template <typename T>
	bool operator()(const T * const parameter_0, T *residual) const
	{
		T ro, pt, he;
		ro = parameter_0[0];
		pt = parameter_0[1];
		he = parameter_0[2];

		T r1, r2, r3, r4, r5, r6, r7, r8, r9;
		r1 = cos(he)*cos(ro) + sin(he)*sin(pt)*sin(ro);
		r2 = cos(he)*sin(pt)*sin(ro) - cos(ro)*sin(he);
		r3 = -cos(pt)*sin(ro);
		r4 = cos(pt)*sin(he);
		r5 = cos(he)*cos(pt);
		r6 = sin(pt);
		r7 = cos(he)*sin(ro) - cos(ro)*sin(he)*sin(pt);
		r8 = -sin(he)*sin(ro) - cos(he)*cos(ro)*sin(pt);
		r9 = cos(pt)*cos(ro);

		T bx, by, bz;
		bx = _x - T(_xs);
		by = _y - T(_ys);
		bz = _z - T(_zs);

		T A, B, C;
		A = r1*bx + r2*by + r3*bz;
		B = r4*bx + r5*by + r6*bz;
		C = r7*bx + r8*by + r9*bz;

		T tanLon = atan2(A, B);
		T tanLat = atan2(sqrt(A * A + B * B), C);

		residual[0] = tanLon - T(_px);
		residual[1] = tanLat - T(_py);

		return true;
	}

private:
	const double _px, _py, _x, _y, _z;
	const double _xs, _ys, _zs;
};

// 只解算roll/pitch/roll
int SPP_S::solvePanoParameter_ceres3(panoPara &pp, std::vector<pointData> pd)
{
	double para_0[3] = { pp.xs, pp.ys, pp.zs };
	double para_1[3] = { pp.alpha, pp.phi, pp.belta };

	ceres::Problem problem;
	for (int i = 0; i < pd.size(); i++)
	{
		pd[i].px = pd[i].px * DPI - PI;
		pd[i].py = pd[i].py * DPI;

		problem.AddResidualBlock(new ceres::AutoDiffCostFunction<ResectionResidual3, 2, 3>
			(new ResectionResidual3(pd[i].px, pd[i].py, pd[i].x, pd[i].y, pd[i].z, pp.xs, pp.ys, pp.zs)),
			NULL, /*para_0,*/ para_1);
	}

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	ceres::Solve(options, &problem, &summary);

	pp.xs = para_0[0];
	pp.ys = para_0[1];
	pp.zs = para_0[2];
	pp.alpha = para_1[0];
	pp.phi = para_1[1];
	pp.belta = para_1[2];

	return 0;
}