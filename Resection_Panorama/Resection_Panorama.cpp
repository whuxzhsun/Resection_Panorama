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

void calib_L10();

int _tmain(int argc, _TCHAR* argv[])
{
	calib_L10();

	return 0;
}

void calib_L10()
{
	vector<panoPara>	pps;

	/*
		读取图像参数
	*/
	string camFile("G:\\20170204_L10\\20180204_BEIJING_L10_TEST_02\\image\\pixel\\imageList_no_offset.txt");
	string pointsFile("G:\\20170204_L10\\20180204_BEIJING_L10_TEST_02\\image\\pixel\\pixel_points.txt");

	ifstream readCam(camFile, ios::in);
	if (!readCam.is_open())
	{
		cout << "failed to open file:" << camFile << endl;
		return;
	}

	while (!readCam.eof())
	{
		char lineStr[512];
		if (!readCam.getline(lineStr, 512))
			break;

		double x, y, z, al, ph, de;
		string second, date, time;
		double la, lo;

		panoPara pp;
		stringstream ss;
		ss << lineStr << endl;
		ss >> pp.imgName >> la >> lo >> x >> y >> z >> al >> ph >> de >> second >> date >> time;
		pp.xs = x;		pp.ys = y;	pp.zs = z;
		pp.alpha = al;	pp.phi = ph;	pp.belta = de;

		pps.push_back(pp);
	}

	readCam.close();

	/*
		读取每幅图像对应的像点
	*/
	ifstream readPoints(pointsFile, ios::in);
	if (!readPoints.is_open())
	{
		cout << "failed to open file:" << pointsFile << endl;
		return;
	}

	while (!readPoints.eof())
	{
		char lineStr[512];
		if (!readPoints.getline(lineStr, 512))
			break;

		pointData pd;
		stringstream ss;
		ss << lineStr << endl;
		ss >> pd.imgName >> pd.px >> pd.py >> pd.x >> pd.y >> pd.z;

		for (int i = 0; i < pps.size(); i++)
		{
			if (pd.imgName == pps[i].imgName)
			{
				pps[i].pixels.push_back(pd);
			}
		}
	}

	readPoints.close();

	SPP_S spp;
	spp.setImgSize(8192, 4096);
	spp.solvePanoParameter(pps);


/*	ifstream inCammer("D:\\Data\\20180203_BEIJING_L10\\image\\resection\\75.txt", ios::in);
	char fileHeader[512];
	inCammer.getline(fileHeader, 512);
	ifstream inPoint("D:\\Data\\20180203_BEIJING_L10\\image\\resection\\ptpx.txt", ios::in);

	ofstream outCam("D:\\Data\\20180203_BEIJING_L10\\image\\resection\\Calibration.txt", ios::out);
	outCam << fileHeader << endl;
	ofstream outContrast("D:\\Data\\20180203_BEIJING_L10\\image\\resection\\Contrast.txt", ios::out);
	outContrast << "old_X\t" << "old_Y\t" << "old_Z\t" << "old_roll\t" << "old_pitch\t" << "old_heanding\t"
		<< "new_X\t" << "new_Y\t" << "new_Z\t" << "new_roll\t" << "new_pitch\t" << "new_heanding\t" << endl;

	ofstream outDiff("D:\\Data\\20180203_BEIJING_L10\\image\\resection\\Diff.txt", ios::out);

	ofstream outMeanError("D:\\Data\\20180203_BEIJING_L10\\image\\resection\\MeanError.txt", ios::out);

	outCam.setf(ios::fixed);
	outCam.width(16);

	outContrast.setf(ios::fixed);
	outContrast.width(16);

	outDiff << "mean_x		mean_y		mean_z		mean_roll	mean_pitch	mean_yaw	diff_x		diff_y		diff_z		diff_roll	diff_pitch	diff_yaw\n";
	outDiff.setf(ios::fixed);
	outDiff.width(16);

	outMeanError.setf(ios::fixed);
	outMeanError.width(16);

	for (k = 0; k < 1; k++)
	{
		panoPara pp;
		vector<pointData> pd;

		int n, px, py;
		double x, y, z, al, ph, de;
		string imgName, second, date, time;
		double la, lo;

		inCammer >> imgName >> la >> lo >> x >> y >> z >> al >> ph >> de >> second >> date >> time;
		pp.xs = x;		pp.ys = y;	pp.zs = z;
		pp.alpha = al;	pp.phi = ph;	pp.belta = de;
		panoPara oldPP = pp;

		pointData pt;
		for (int i = 0; i < 5; i++)
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