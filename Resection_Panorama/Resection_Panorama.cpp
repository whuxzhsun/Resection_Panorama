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

int getRotation(double inRoll, double inPitch, double inYaw, double outR[])
{
	inRoll = -inRoll;
	inPitch = -inPitch;
	outR[0] = cos(inRoll)*cos(inYaw) + sin(inPitch)*sin(inRoll)*sin(inYaw);
	outR[1] = cos(inPitch)*sin(inYaw);
	outR[2] = -sin(inRoll)*cos(inYaw) + sin(inPitch)*cos(inRoll)*sin(inYaw);
	outR[3] = -cos(inRoll)*sin(inYaw) + sin(inPitch)*sin(inRoll)*cos(inYaw);
	outR[4] = cos(inPitch)*cos(inYaw);
	outR[5] = sin(inRoll)*sin(inYaw) + sin(inPitch)*cos(inRoll)*cos(inYaw);
	outR[6] = cos(inPitch)*sin(inRoll);
	outR[7] = -sin(inPitch); // + or - ?
	outR[8] = cos(inRoll)*cos(inPitch);

	return 0;
}

void resection_L10()
{
	ifstream inCammer("pts\\img50.txt", ios::in);
	char fileHeader[512];
//	inCammer.getline(fileHeader, 512);
	ifstream inPoint("pts\\50.txt", ios::in);
	int nPts = 4;

	ofstream outCam("pts\\Calibration.txt", ios::out);
	outCam << fileHeader << endl;
	ofstream outContrast("pts\\Contrast.txt", ios::out);
	outContrast << "old_X\t" << "old_Y\t" << "old_Z\t" << "old_roll\t" << "old_pitch\t" << "old_heanding\t"
		<< "new_X\t" << "new_Y\t" << "new_Z\t" << "new_roll\t" << "new_pitch\t" << "new_heanding\t" << endl;

	ofstream outDiff("pts\\Diff.txt", ios::out);

	ofstream outMeanError("pts\\MeanError.txt", ios::out);

	outCam.setf(ios::fixed);
	outCam.width(16);

	outContrast.setf(ios::fixed);
	outContrast.width(16);

	outDiff << "mean_x		mean_y		mean_z		mean_roll	mean_pitch	mean_yaw	diff_x		diff_y		diff_z		diff_roll	diff_pitch	diff_yaw\n";
	outDiff.setf(ios::fixed);
	outDiff.width(16);

	outMeanError.setf(ios::fixed);
	outMeanError.width(16);

	for (int k = 0; k < 1; k++)
	{
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
			inPoint >> n >> px >> py  >> x >> y >> z;
			pt.px = px;	pt.py = py;
			pt.x = x;	pt.y = y;	pt.z = z;
			pd.push_back(pt);
		}

		SPP_S spp(8192, 4096);
		spp.solvePanoParameter(pp, pd);

		outCam << imgName << "\t" << la << "\t" << lo << "\t" << pp.xs << "\t" << pp.ys << "\t"
			<< pp.zs << "\t" << pp.alpha << "\t" << pp.phi << "\t" << pp.belta << "\t"
			<< second << "\t" << date << "\t" << time << endl;

		outContrast << oldPP.xs << "\t" << oldPP.ys << "\t" << oldPP.zs << "\t" << oldPP.alpha << "\t" 
			<< oldPP.phi << "\t" << oldPP.belta << "\t" << pp.xs << "\t" << pp.ys << "\t"
			<< pp.zs << "\t" << pp.alpha << "\t" << pp.phi << "\t" << pp.belta << endl;

		outDiff << pp.mean0[0] << "\t" << pp.mean0[1] << "\t" <<
			pp.mean0[2] << "\t" <<	pp.mean0[3] << "\t" <<
			pp.mean0[4] << "\t" <<	pp.mean0[5] << "\t" 
			<< oldPP.xs - pp.xs << "\t" << oldPP.ys - pp.ys << "\t" << oldPP.zs - pp.zs << "\t"
			<< oldPP.alpha - pp.alpha << "\t" << oldPP.phi - pp.phi << "\t" << oldPP.belta - pp.belta << "\n";

		outMeanError << pp.xs << "\t" << pp.mean0[0] << "\n" <<
			pp.ys	<< "\t" << pp.mean0[1] << "\n" <<
			pp.zs	<< "\t" << pp.mean0[2] << "\n" <<
			pp.alpha<< "\t" << pp.mean0[3] << "\n" <<
			pp.phi	<< "\t" << pp.mean0[4] << "\n" <<
			pp.belta<< "\t" << pp.mean0[5] << "\n\n" ;

		double dd[3] = { 0.0864, -0.2961, 0.3015 };

		double R[9];
		getRotation(pp.alpha/* * DEG_TO_RAD*/, pp.phi/* * DEG_TO_RAD*/, pp.belta/* * DEG_TO_RAD*/, R);

		double ax = R[0] * dd[0] + R[1] * dd[1] + R[2] * dd[2];
		double ay = R[3] * dd[0] + R[4] * dd[1] + R[5] * dd[2];
		double az = R[6] * dd[0] + R[7] * dd[1] + R[8] * dd[2];

		outMeanError << 01 << " " << second + 24 * 3600 << "  " << pp.xs - ax << "  " << pp.ys - ay << "  " << pp.zs - az << "  "
			<< "  0" << "  0" << "  0" << "  0" << "  0" << "  0\n";

		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.xs << "\t" ;
		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.mean0[0] << "\n";
		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.ys << "\t" ;
		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.mean0[1] << "\n";
		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.zs << "\t" ;
		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.mean0[2] << "\n";
		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.alpha << "\t" ;
		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.mean0[3] << "\n";
		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.phi << "\t" ;
		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.mean0[4] << "\n";
		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.belta << "\t" ;
		cout.width(12);
		cout.setf(ios::fixed);
		cout.precision(4);
		cout << pp.mean0[5] << "\n";

		cout<< oldPP.xs - pp.xs << "\t" << oldPP.ys - pp.ys << "\t" << oldPP.zs - pp.zs << "\n"
			<< oldPP.alpha - pp.alpha << "\t" << oldPP.phi - pp.phi << "\t" << oldPP.belta - pp.belta << "\n\n";

	}
	inCammer.close();
	outCam.close();
	outContrast.close();
	inPoint.close();
	outDiff.close();
	outMeanError.close();//*/
}