/*先頭で日本語を打ち込んでおけばソースツリーで表示したときに文字化けしないらしいので*/
#ifndef __TEMPLATE_MATCHONG__
#define __TEMPLATE_MATCHONG__

#include "narivectorpp.h"
//#include "nmi_matching_It.h"
#include <omp.h>
#include "stdlib.h"

template <class T, class M>
void template_mathcing(const nari::vector<M> &imgRef, const nari::vector<T> &imgFl, nari::vector<nari::vector<int>> &DispRef,
	const nari::vector<nari::vector<int>> &DispFl, int xeRef, int yeRef, int zeRef, int xeFl, int yeFl, int zeFl,
	int tmp, int rangex, int rangey, int rangez)
{
	int disp_num = DispFl.size();
	std::cout << disp_num << std::endl;

	int tmp_size = (tmp * 2 + 1)*(tmp * 2 + 1)*(tmp * 2 + 1);

	//読み見込んだ計測点すべてを走査するループ
	//※※NMIの関数を使うとなぜか並列化で値が混ざるのでNMIを使用する際は並列化は外してください※※
//#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int a = 0; a < disp_num; a++) {
		nari::vector<int> disp(3);
		nari::vector<short> tmp_Fl(tmp_size);
		nari::vector<int> tmp_space((tmp_size), 0);
		int t = 0;
		int v = 0;
		//テンプレートマッチング前にあらかじめ位置合わせしたい側の画像のテンプレート作成
		for (int r = 0; r < 2 * tmp + 1; r++) {
			for (int q = 0; q < 2 * tmp + 1; q++) {
				for (int p = 0; p < 2 * tmp + 1; p++) {
					int x = DispFl[a][0] - tmp + p;
					int y = DispFl[a][1] - tmp + q;
					int z = DispFl[a][2] - tmp + r;
					//テンプレートが画像からはみ出た場合(0,0,0)の濃度値を入れる
					if ((x >= 0) && (y >= 0) && (z >= 0) && (x < xeFl) && (y < yeFl) && (z < zeFl)) {
						int s = xeFl*yeFl*z + xeFl*y + x;
						//テンプレート内の画素を0〜32に正規化
						//tmp_Ref[t] = imgRef[s]*32/65535;
						tmp_Fl[t] = imgFl[s];
						t++;
					}
					else {
						//tmp_spaceにテンプレートの画素が入ってない画素番号を１に
						tmp_Fl[t] = 500;
						tmp_space[t] = 1;
						t++;
					}
				}
			}
		}

		std::cout << "(^^)<テンプレート作った" << std::endl;


		//相関係数を格納する変数を定義
		double cc = 0, cc_max = 0;

		//対応点の座標を入れる変数
		int xs, ys, zs;

		//テンプレートマッチング開始
		for (int k = 0; k < rangez * 2 + 1; k++) {
			for (int j = 0; j < rangey * 2 + 1; j++) {
				for (int i = 0; i < rangex * 2 + 1; i++) {
					nari::vector<short> tmp_Ref(tmp_size);
					int u = 0;
					//位置合わせされる側の画像テンプレート作成
					for (int r = 0; r < 2 * tmp + 1; r++) {
						for (int q = 0; q < 2 * tmp + 1; q++) {
							for (int p = 0; p < 2 * tmp + 1; p++) {
								int x = DispFl[a][0] - rangex + i - tmp + p;
								int y = DispFl[a][1] - rangey + j - tmp + q;
								int z = DispFl[a][2] - rangez + k - tmp + r;
								//テンプレートが画像からはみ出なければテンプレートつくる
								if ((x >= 0) && (y >= 0) && (z >= 0) && (x < xeRef) && (y < yeRef) && (z < zeRef)) {
									int s = xeRef*yeRef*z + xeRef*y + x;
									//tmp_Fl[u] = imgFl[s]*32/65535;
									tmp_Ref[u] = imgRef[s];
									u++;
								}
								else {
									tmp_Ref[u] = 500;
									u++;
								}
							}
						}
					}

					//ここからテンプレート同士の相関係数を計算
					double meanref = 0.0, meanfl = 0.0;
					int tmp_true = 0;
					for (int c = 0; c < tmp_size; c++) {
						if (tmp_space[c] == 0) {
							meanref += tmp_Ref[c];
							meanfl += tmp_Fl[c];
							tmp_true++;
						}
					}
					meanref = meanref / tmp_true;
					meanfl = meanfl / tmp_true;

					double stdref = 0.0, stdfl = 0.0, cov = 0.0;
					for (int c = 0; c < tmp_size; c++) {
						if (tmp_space[c] == 0) {
							stdref += (tmp_Ref[c] - meanref)*(tmp_Ref[c] - meanref);
							stdfl += (tmp_Fl[c] - meanfl)*(tmp_Fl[c] - meanfl);
							cov += (tmp_Ref[c] - meanref)*(tmp_Fl[c] - meanfl);
						}
					}

					//分散0になると分母が0になるのでスキップ
					if (stdfl != 0) {
						cc = cov / (sqrt(stdref)*sqrt(stdfl));
					}

					//相関係数が最大値をとるときの座標を保存

					if ((cc >= cc_max) && (DispFl[a][0] - rangex + i >= 0) && (DispFl[a][1] - rangey + j >= 0) && (DispFl[a][2] - rangez + k >= 0)) {
						cc_max = cc;
						xs = DispFl[a][0] - rangex + i;
						ys = DispFl[a][1] - rangey + j;
						zs = DispFl[a][2] - rangez + k;
					}
				}
			}
		}
		disp[0] = xs;
		disp[1] = ys;
		disp[2] = zs;
		std::cout << "x=" << DispFl[a][0] << ", y=" << DispFl[a][1] << ", z=" << DispFl[a][2] << std::endl;
		std::cout << "xs=" << xs << ", ys=" << ys << ", zs=" << zs << std::endl;
		DispRef.push_back(disp);

	}
}
#endif