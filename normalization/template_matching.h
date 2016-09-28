/*�擪�œ��{���ł�����ł����΃\�[�X�c���[�ŕ\�������Ƃ��ɕ����������Ȃ��炵���̂�*/
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

	//�ǂ݌����񂾌v���_���ׂĂ𑖍����郋�[�v
	//����NMI�̊֐����g���ƂȂ������񉻂Œl��������̂�NMI���g�p����ۂ͕��񉻂͊O���Ă�����������
//#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int a = 0; a < disp_num; a++) {
		nari::vector<int> disp(3);
		nari::vector<short> tmp_Fl(tmp_size);
		nari::vector<int> tmp_space((tmp_size), 0);
		int t = 0;
		int v = 0;
		//�e���v���[�g�}�b�`���O�O�ɂ��炩���߈ʒu���킹���������̉摜�̃e���v���[�g�쐬
		for (int r = 0; r < 2 * tmp + 1; r++) {
			for (int q = 0; q < 2 * tmp + 1; q++) {
				for (int p = 0; p < 2 * tmp + 1; p++) {
					int x = DispFl[a][0] - tmp + p;
					int y = DispFl[a][1] - tmp + q;
					int z = DispFl[a][2] - tmp + r;
					//�e���v���[�g���摜����͂ݏo���ꍇ(0,0,0)�̔Z�x�l������
					if ((x >= 0) && (y >= 0) && (z >= 0) && (x < xeFl) && (y < yeFl) && (z < zeFl)) {
						int s = xeFl*yeFl*z + xeFl*y + x;
						//�e���v���[�g���̉�f��0�`32�ɐ��K��
						//tmp_Ref[t] = imgRef[s]*32/65535;
						tmp_Fl[t] = imgFl[s];
						t++;
					}
					else {
						//tmp_space�Ƀe���v���[�g�̉�f�������ĂȂ���f�ԍ����P��
						tmp_Fl[t] = 500;
						tmp_space[t] = 1;
						t++;
					}
				}
			}
		}

		std::cout << "(^^)<�e���v���[�g�����" << std::endl;


		//���֌W�����i�[����ϐ����`
		double cc = 0, cc_max = 0;

		//�Ή��_�̍��W������ϐ�
		int xs, ys, zs;

		//�e���v���[�g�}�b�`���O�J�n
		for (int k = 0; k < rangez * 2 + 1; k++) {
			for (int j = 0; j < rangey * 2 + 1; j++) {
				for (int i = 0; i < rangex * 2 + 1; i++) {
					nari::vector<short> tmp_Ref(tmp_size);
					int u = 0;
					//�ʒu���킹����鑤�̉摜�e���v���[�g�쐬
					for (int r = 0; r < 2 * tmp + 1; r++) {
						for (int q = 0; q < 2 * tmp + 1; q++) {
							for (int p = 0; p < 2 * tmp + 1; p++) {
								int x = DispFl[a][0] - rangex + i - tmp + p;
								int y = DispFl[a][1] - rangey + j - tmp + q;
								int z = DispFl[a][2] - rangez + k - tmp + r;
								//�e���v���[�g���摜����͂ݏo�Ȃ���΃e���v���[�g����
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

					//��������e���v���[�g���m�̑��֌W�����v�Z
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

					//���U0�ɂȂ�ƕ��ꂪ0�ɂȂ�̂ŃX�L�b�v
					if (stdfl != 0) {
						cc = cov / (sqrt(stdref)*sqrt(stdfl));
					}

					//���֌W�����ő�l���Ƃ�Ƃ��̍��W��ۑ�

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