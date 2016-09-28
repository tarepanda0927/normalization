/*先頭で日本語を打ち込んでおけばソースツリーで表示したときに文字化けしないらしいので*/
#include "naricommon.h"
#include "ip/narigaussian.h"
#include <iomanip>

#include "naritimer.h"
#include "naricaseinfo.h"
#include <algorithm>
#include "other/narimhd.h"

#include <sstream>
#include "naricommon.h"
#include "narisystem.h"
#include "naripathline.h"
#include "naricaseinfo.h"
#include "ip/narimorphology.h"
#include "ip/narilabeling.h"

#include "ip/naricontour.h"
#include "ip/naridistance.h"
#include "other/narimhd.h"
#include "naritimer.h"
#include "ip/narirbf.h"
#include "../mist/vector.h"
#include "narivectorpp.h"
#include "ip/nariinterpolate.h"
#include "info.h"
#include "raw_io.h"
#include "template_matching.h"
#include "set_point.h"

#include <chrono>
#include <cmath>
#include <iostream>

template <class IMG_T>
void Img_deformation_using_movement(int xeRef, int yeRef, int zeRef, int xeFl, int yeFl, int zeFl, int move_x, int move_y, int move_z, nari::vector<float> &imgMoveX, nari::vector<float> &imgMoveY, nari::vector<float> &imgMoveZ, nari::vector<IMG_T> &imgI, nari::vector<IMG_T> &imgO)
{

	for (int z = 0; z < zeRef; z++)
	{
		for (int y = 0; y < yeRef; y++)
		{
			for (int x = 0; x < xeRef; x++)
			{
				int s = z * xeRef * yeRef + y * xeRef + x;
				//フチに500が入るように設定
				if ((-1 < (x - move_x)) && ((x - move_x) < xeFl) && (-1 < (y - move_y)) && ((y - move_y) < yeFl) && (-1 < (z - move_z)) && ((z - move_z) < zeFl)) {
					imgO[s] = static_cast<IMG_T>(nari::interpolate_value::linear(imgI.ptr(), imgMoveX[s], imgMoveY[s], imgMoveZ[s], xeFl, yeFl, zeFl));
				}
				else {

					imgO[s] = 500;
				}
			}
		}
	}
}
template <class IMG_T>
void Label_deformation_using_movement(int xeRef, int yeRef, int zeRef, int xeFl, int yeFl, int zeFl, nari::vector<float> &imgMoveX, nari::vector<float> &imgMoveY, nari::vector<float> &imgMoveZ, nari::vector<IMG_T> &imgI, nari::vector<IMG_T> &imgO)
{
	for (int z = 0; z < zeRef; z++)
	{
		for (int y = 0; y < yeRef; y++)
		{
			for (int x = 0; x < xeRef; x++)
			{
				int s = z * xeRef * yeRef + y * xeRef + x;
				imgO[s] = static_cast<IMG_T>(nari::interpolate_value::linear(imgI.ptr(), imgMoveX[s], imgMoveY[s], imgMoveZ[s], xeFl, yeFl, zeFl));
			}
		}
	}
}

template <typename T>
T median(std::vector<T>& c)
{
	if (c.size() % 2 == 0) {
		size_t n1 = floor(c.size() / 2);
		size_t n2 = round(c.size() / 2);

	}
	else {
		size_t n = c.size() / 2;
		std::nth_element(c.begin(), c.begin() + n, c.end());
		return c[n];
	}
}
void main(int argc, char *argv[])
{

	info input_info;
	input_info.input(argv[1]);

	//テキストデータ読み込み
	std::vector<std::string> rcase;
	std::ifstream r_case(input_info.dir_list + input_info.case_rlist);
	std::string buf_ft;
	std::string buf_rt;
	while (r_case&& getline(r_case, buf_rt))
	{
		rcase.push_back(buf_rt);
	}

	for (int i = 0; i < rcase.size(); i++) {

		//テスト症例読み込み
		nari::mhd mhdI, mhdIL;
		mhdI.load(input_info.dir_Ref + rcase[i] + ".mhd");
		mhdIL.load(input_info.dir_Ref + rcase[i] + "_label.mhd");
		int xeRef = mhdI.size1();
		int yeRef = mhdI.size2();
		int zeRef = mhdI.size3();
		double xrRef = mhdI.reso1();
		double yrRef = mhdI.reso2();
		double zrRef = mhdI.reso3();
		nari::vector<short> imgI;
		nari::vector<unsigned char> imgIL;
		imgI.load_file_bin(input_info.dir_Ref + rcase[i] + ".raw");
		imgIL.load_file_bin(input_info.dir_Ref + rcase[i] + "_label.raw");

		//テンプレートマッチングパラメータ読み込み
		int rx = input_info.rangex;
		int ry = input_info.rangey;
		int rz = input_info.rangez;
		int tmp = input_info.tmp;
		int tmp_size = (tmp * 2 + 1)*(tmp * 2 + 1)*(tmp * 2 + 1);
		//std::vector<double> a_x, a_y, a_z, b_x, b_y, b_z;
		double a_x = 0, a_y = 0, a_z = 0, b_x = 0, b_y = 0, b_z = 0;
		for (int j = 0; j < rcase.size(); j++) {
			if (j == i)continue;  //テスト症例はスキップ
			//正解ラベルを使ってきれいに位置合わせしてある学習症例を読み込む
			nari::vector<short> imgA;
			nari::vector<unsigned char> imgAL;
			imgA.load_file_bin(input_info.dir_mvd + rcase[j] + ".raw");
			imgAL.load_file_bin(input_info.dir_Ans + rcase[j] + "_label.raw");
			//浮動症例のラベル表面にランドマーク自動決定
			std::cout << "setpoint" << std::endl;
			nari::vector<nari::vector<int>> DispA;
			nari::vector<nari::vector<int>> preA;
			std::string dir_Alist = input_info.dir_out + "setpoint.txt";
			set_point(dir_Alist, imgAL, xeRef, yeRef, zeRef, DispA);

			//[26][27]の中点と[0]の座標をpreFlに追加
			int x_46 = (DispA[26][0] + DispA[27][0]) / 2;
			int y_46 = DispA[26][1];
			int z_46 = DispA[26][2];
			int x_0 = DispA[0][0];
			int y_0 = DispA[0][1];
			int z_0 = DispA[0][2];
			std::cout << DispA[0][0] << "," << DispA[0][1] << "," << DispA[0][2] << std::endl;
			//preFlに保存
			nari::vector<int> disp(3);
			disp[0] = x_46;
			disp[1] = y_46;
			disp[2] = z_46;
			preA.push_back(disp);
			disp[0] = x_0;
			disp[1] = y_0;
			disp[2] = z_0;
			preA.push_back(disp);

			//std::cout << preA[0][0] << " ." << preA[0][1] << " ." << preA[0][2] << std::endl;
			//std::cout << preA[1][0] << " ." << preA[1][1] << " ." << preAl[1][2] << std::endl;

			////ここからの処理は大まかな平行移動を行う処理
			nari::vector<nari::vector<int>> preI;
			//2点のみテンプレートマッチング
			template_mathcing(imgI, imgA, preI, preA, xeRef, yeRef, zeRef,
				xeRef, yeRef, zeRef, tmp, rx, ry, rz);
			a_x += preI[0][0];
			a_y += preI[0][1];
			a_z += preI[0][2];
			b_x += preI[1][0];
			b_y += preI[1][1];
			b_z += preI[1][2];
			/*a_x.push_back(preI[0][0]);
			a_y.push_back(preI[0][1]);
			a_z.push_back(preI[0][2]);
			b_x.push_back(preI[1][0]);
			b_y.push_back(preI[1][1]);
			b_z.push_back(preI[1][2]);*/
		}
		
		nari::vector<nari::vector<int>> preI;
		nari::vector<int> disp(3);
		disp[0] = round(a_x / (rcase.size() - 1));
		disp[1] = round(a_y / (rcase.size() - 1));
		disp[2] = round(a_z / (rcase.size() - 1));
		preI.push_back(disp);
		disp[0] = round(b_x / (rcase.size() - 1));
		disp[1] = round(b_y / (rcase.size() - 1));
		disp[2] = round(b_z / (rcase.size() - 1));
		preI.push_back(disp);
		/*std::sort(a_x.begin(), a_x.end());
		std::sort(a_y.begin(), a_y.end());
		std::sort(a_z.begin(), a_z.end());
		std::sort(b_x.begin(), b_x.end());
		std::sort(b_y.begin(), b_y.end());
		std::sort(b_z.begin(), b_z.end());*/
		/*int n = a_x.size();
		if (n % 2 == 0) {
			disp[0] = round((a_x[n / 2] + a_x[n / 2 + 1]) / 2);
			disp[1] = round((a_y[n / 2] + a_y[n / 2 + 1]) / 2);
			disp[2] = round((a_z[n / 2] + a_z[n / 2 + 1]) / 2);
			preI.push_back(disp);
			disp[0] = round((b_x[n / 2] + b_x[n / 2 + 1]) / 2);
			disp[1] = round((b_y[n / 2] + b_y[n / 2 + 1]) / 2);
			disp[2] = round((b_z[n / 2] + b_z[n / 2 + 1]) / 2);
			preI.push_back(disp);
		}
		else {
			int m = n / 2 + 1;
			disp[0] = (int)a_x[m];
			disp[1] = (int)a_y[m];
			disp[2] = (int)a_z[m];
			preI.push_back(disp);
			disp[0] = (int)b_x[m];
			disp[1] = (int)b_y[m];
			disp[2] = (int)b_z[m];
			preI.push_back(disp);
		}*/
		
		std::cout <<"/////result/////" << std::endl;
		std::cout << preI[0][0] << "," << preI[0][1] << "," << preI[0][2] << std::endl;
		std::cout << preI[1][0] << "," << preI[1][1] << "," << preI[1][2] << std::endl;
		nari::vector<float> imgMoveX(xeRef * yeRef * zeRef);
		nari::vector<float> imgMoveY(xeRef * yeRef * zeRef);
		nari::vector<float> imgMoveZ(xeRef * yeRef * zeRef);
		int move_x = 96 - preI[0][0];
		int move_y = 96 - preI[0][1];
		int move_z = 45 - preI[0][2];


		nari::vector<nari::vector<int>> Disp_pre_Fl;


		for (int z = 0; z < zeRef; z++)
		{
			for (int y = 0; y < yeRef; y++)
			{
				for (int x = 0; x < xeRef; x++)
				{
					int s = z * xeRef * yeRef + y * xeRef + x;
					imgMoveX[s] = (float)x - move_x;
					imgMoveY[s] = (float)y - move_y;
					imgMoveZ[s] = (float)z - move_z;
				}
			}
		}

		//簡単な平行移動を行いった結果を保存
		nari::vector<short> imgI2(xeRef * yeRef *zeRef);
		nari::vector<unsigned char>imgIL2(xeRef * yeRef *zeRef);
		Img_deformation_using_movement(xeRef, yeRef, zeRef, xeRef, yeRef, zeRef, move_x, move_y, move_z, imgMoveX, imgMoveY, imgMoveZ, imgI, imgI2);
		Label_deformation_using_movement(xeRef, yeRef, zeRef, xeRef, yeRef, zeRef, imgMoveX, imgMoveY, imgMoveZ, imgIL, imgIL2);

		std::cout << "回転" << std::endl;
		//回転して向きを合わせる
		float refy, refz, cosR, sinR, normR;
		refy = (int)(preI[1][1] + move_y - 96);
		refz = (int)(preI[1][2] + move_z - 45);
		normR = (double)sqrt(refy*refy + refz*refz);
		cosR = (float)(preI[1][1] + move_y - 96) / normR;
		sinR = (float)(preI[1][2] + move_z - 45) / normR;

		for (int z = 0; z < zeRef; z++)
		{
			for (int y = 0; y < yeRef; y++)
			{
				for (int x = 0; x < xeRef; x++)
				{
					int s = z * xeRef * yeRef + y * xeRef + x;
					imgMoveX[s] = (float)x;
					imgMoveY[s] = (float)(y - 96)*(-cosR) + (z - 45)*(sinR)+96;
					imgMoveZ[s] = (float)(y - 96)*(-sinR) + (z - 45)*(-cosR) + 45;
				}
			}
		}

		nari::vector<short> imgI3(xeRef * yeRef * zeRef);
		nari::vector<unsigned char> imgIL3(xeRef * yeRef * zeRef);
		//回転移動を行いった結果をimgFl2に保存
		Label_deformation_using_movement(xeRef, yeRef, zeRef, xeRef, yeRef, zeRef, imgMoveX, imgMoveY, imgMoveZ, imgIL2, imgIL3);
		Label_deformation_using_movement(xeRef, yeRef, zeRef, xeRef, yeRef, zeRef, imgMoveX, imgMoveY, imgMoveZ, imgI2, imgI3);
		mhdI.save_mhd_and_image(imgI3, input_info.dir_out + rcase[i] + ".raw");
		mhdIL.save_mhd_and_image(imgIL3, input_info.dir_out + rcase[i] + "_label.raw");
	}
}