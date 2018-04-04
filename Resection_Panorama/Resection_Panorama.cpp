// Resection_Panorama.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "solvePanoPara.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
/*#include <iomanip>*/

using namespace std;

#define DEG_TO_RAD	.017453292519943296

void resection_L10();

int _tmain(int argc, _TCHAR* argv[])
{
	resection_L10();

	return 0;
}

void resection_L10()
{
	int nPts = 15;

	ifstream inCammer("pts\\img50.txt", ios::in);
	ifstream inPoint("pts\\50.txt", ios::in);
	ofstream outCam("pts\\Calibration.txt", ios::out);
	ofstream outDiff("pts\\Diff.txt", ios::out);
	ofstream outMeanError("pts\\MeanError.txt", ios::out);

	outCam.setf(ios::fixed);
	outCam.width(16);

//	outDiff << "diff_x		diff_y		diff_z		diff_roll	diff_pitch	diff_yaw\n";
	outDiff.setf(ios::fixed);
	outDiff.width(16);

	outMeanError.setf(ios::fixed);
	outMeanError.width(16);

	panoPara pp;
	vector<pointData> pd;

	int n, px, py;
	double x, y, z, al, ph, de;
	string imgName, date, time;
	double la, lo, second;

	inCammer >> imgName >> la >> lo >> x >> y >> z >> al >> ph >> de >> second >> date >> time;
	pp.xs = x;		pp.ys = y;	pp.zs = z;
	pp.alpha = al;	pp.phi = ph;	pp.belta = de;
	panoPara oldPP = pp;

	pointData pt;
	for (int i = 0; i < nPts; i++)
	{
		inPoint >> n >> px >> py >> x >> y >> z;
		pt.px = px;	pt.py = py;
		pt.x = x;	pt.y = y;	pt.z = z;
		pd.push_back(pt);
	}

	panoPara truePP;
	truePP.xs = 458317.6695;
	truePP.ys = 4403414.6670;
	truePP.zs = 23.3259;
	truePP.alpha = 0.015362;
	truePP.phi   = 0.053639;
	truePP.belta = 3.173437;

	double pixelToAngle = 360.0 / 8192;

	pointData pt01, pt02, pt03;
	for (int i = 0; i < nPts; i++)
	{
		pt01 = pd[i];
		for (int j = i + 1; j < nPts; j++)
		{
			pt02 = pd[j];
			for (int k = j + 1; k < nPts; k++)
			{
				pt03 = pd[k];

				vector<pointData> pd2;
				pd2.push_back(pt01);
				pd2.push_back(pt02);
				pd2.push_back(pt03);

				pp = oldPP;

				SPP_S spp(8192, 4096);
				spp.solvePanoParameter(pp, pd2);

				outCam << imgName << "\t" << la << "\t" << lo << "\t" << pp.xs << "\t" << pp.ys << "\t"
					<< pp.zs << "\t" << pp.alpha << "\t" << pp.phi << "\t" << pp.belta << "\t"
					<< second << "\t" << date << "\t" << time << endl;

// 				outDiff << oldPP.xs - pp.xs << "\t" << oldPP.ys - pp.ys << "\t" << oldPP.zs - pp.zs << "\t"
// 					<< oldPP.alpha - pp.alpha << "\t" << oldPP.phi - pp.phi << "\t" << oldPP.belta - pp.belta << "\n";

				outMeanError << pp.xs << "\t" << pp.mean0[0] << "\n" <<
					pp.ys << "\t" << pp.mean0[1] << "\n" <<
					pp.zs << "\t" << pp.mean0[2] << "\n" <<
					pp.alpha << "\t" << pp.mean0[3] << "\n" <<
					pp.phi << "\t" << pp.mean0[4] << "\n" <<
					pp.belta << "\t" << pp.mean0[5] << "\n\n";

				outDiff << i << j << k << "\t";
				outDiff.setf(ios::fixed);
				outDiff.width(12);
				outDiff.precision(4);
				outDiff << pixelToAngle * (abs(pt01.px - pt02.px) + abs(pt01.py - pt02.py)) / 1.414 << "\t";
				outDiff.setf(ios::fixed);
				outDiff.width(12);
				outDiff.precision(4);
				outDiff << pixelToAngle * (abs(pt02.px - pt03.px) + abs(pt02.py - pt03.py)) / 1.414 << "\t";
				outDiff.setf(ios::fixed);
				outDiff.width(12);
				outDiff.precision(4);
				outDiff << pixelToAngle * (abs(pt01.px - pt03.px) + abs(pt01.py - pt03.py)) / 1.414 << "\t";
				outDiff.setf(ios::fixed);
				outDiff.width(12);
				outDiff.precision(4);
				outDiff << pp.xs - truePP.xs << "\t";
				outDiff.width(12);
				outDiff.precision(4);
				outDiff << pp.alpha - truePP.alpha << "\t";
				outDiff.width(12);
				outDiff.precision(4);
				outDiff << pp.ys - truePP.ys << "\t";
				outDiff.width(12);
				outDiff.precision(4);
				outDiff << pp.phi - truePP.phi << "\t";
				outDiff.width(12);
				outDiff.precision(4);
				outDiff << pp.zs - truePP.zs << "\t";
				outDiff.width(12);
				outDiff.precision(4);
				outDiff << pp.belta - truePP.belta << "\n";
			}
		}
	}

	/*
	cout.setf(ios::fixed);
	cout.width(12);
	cout.precision(4);
	cout << pp.xs << "\t";
	cout.width(12);
	cout.precision(4);
	cout << pp.mean0[0] << "\n";
	cout.width(12);
	cout.precision(4);
	cout << pp.ys << "\t";
	cout.width(12);
	cout.precision(4);
	cout << pp.mean0[1] << "\n";
	cout.width(12);
	cout.precision(4);
	cout << pp.zs << "\t";
	cout.width(12);
	cout.precision(4);
	cout << pp.mean0[2] << "\n";
	cout.width(12);
	cout.precision(6);
	cout << pp.alpha << "\t";
	cout.width(12);
	cout.precision(4);
	cout << pp.mean0[3] << "\n";
	cout.width(12);
	cout.precision(6);
	cout << pp.phi << "\t";
	cout.width(12);
	cout.precision(4);
	cout << pp.mean0[4] << "\n";
	cout.width(12);
	cout.precision(6);
	cout << pp.belta << "\t";
	cout.width(12);
	cout.precision(4);
	cout << pp.mean0[5] << "\n";

	cout << oldPP.xs - pp.xs << "\t" << oldPP.ys - pp.ys << "\t" << oldPP.zs - pp.zs << "\n"
		<< oldPP.alpha - pp.alpha << "\t" << oldPP.phi - pp.phi << "\t" << oldPP.belta - pp.belta << "\n\n";
		*/

	inCammer.close();
	outCam.close();
	inPoint.close();
	outDiff.close();
	outMeanError.close();
}