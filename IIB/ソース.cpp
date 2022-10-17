#include <iostream>
#include <vector>
#include <intrin.h>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>  

using namespace std;

//�萔�J�n
vector<vector<long long int>> G;//���͂����O���t
long long int n;//���_��
long long int m;//�}��
long long int k;//�p�����[�^(�ŏI�I�Ɋ������钸�_�̐��̍ő�l)
long long int l;//�p�����[�^(�Ɖu���������钸�_�̐��̍ő�l)
vector<long long int> T;//�֐�
long long int x = 5;//�n�b�V���l�̌v�Z�ɗp����hash-multiplier
long long int mod;//�n�b�V���l�̌v�Z�ɗp����mod(���_���ɂ���ĕω�����)
vector<long long int> Monomials;//�n�b�V���l�̌v�Z�ɗp����P����(Monomials[i]��x^i��\��)
//�萔�I��
//�ŏI�I�Ɋ������钸�_�����߂�֐�
long long int who_is_influenced(long long int bit) {
	long long int influenced = 0;//�������Ă��钸�_��bit1�ŕ\��
	for (long long int u = 0; u < n; u++) {
		if (!(bit & (static_cast<long long>(1) << u)) && T[u] == 0)influenced |= (static_cast<long long>(1) << u);
	}
	long long int before = -1;
	while (influenced != before) {
		before = influenced;
		for (long long int u = 0; u < n; u++) {
			if (bit & (static_cast<long long>(1) << u))continue;
			long long int count = 0;
			for (long long unsigned i = 0; i < G[u].size(); i++) {
				if (influenced & (static_cast<long long>(1) << G[u][i]))count++;
			}
			if (count >= T[u])influenced |= (static_cast<long long>(1) << u);
		}
	}
	return influenced;
}

//Monomials�����߂�
void calculate_Monomials() {
	long long int count = 1;
	for (long long int i = 0; i < n; i++) {
		Monomials.push_back(count);
		count = (count*x)%mod;
	}
}
//polynominal hash function(���ꂼ��̒��_�̋ߖT�ɑ΂���n�b�V���l�����߂�)


//���C���֐�

int main() {
	//���͊J�n
	ifstream ifs("paper_sample.txt");

	if (!ifs) {
		std::cout << "Error!";
		return 1;
	}

	string s;
	long long int count = 0;
	long long int x = 0;
	long long int y = 0;

	while (getline(ifs, s, ' ')) {     // �X�y�[�X�i' '�j�ŋ�؂��āC�i�[
		cout << "count:" << count << endl;
		if (count == 0) {//���_��
			n = stoll(s);
			G.resize(n);//�O���t�̑傫���m��
			count++;
		}
		else if (count == 1) {//�}��
			m = stoll(s);
			count++;
		}
		else if (count == 2) {//�p�����[�^k
			k = stoll(s);
			count++;
		}
		else if (count == 3) {//�p�����[�^l
			l = stoll(s);
			count++;
		}
		else if (count > 3 && count < 4 + 2 * m && count % 2 == 0) {//�O���t�̎}�̒[�_
			x = stoll(s);
			x--;
			count++;
		}
		else if (count > 3 && count < 4 + 2 * m && count % 2 == 1) {//�O���t�̎}�̂�����̒[�_
			y = stoll(s);
			y--;
			count++;
			G[x].push_back(y);
			G[y].push_back(x);
		}
		else {//���_��臒l
			x = stoll(s);
			T.push_back(x);
			count++;
		}
	}
	mod = n * n;//mod�͒萔�~n^2
	//���͏I��

	/*���͊m�F�J�n
	cout << "���_��:" << n << endl;
	cout << "�}��:" << m << endl;
	cout << "k:" << k << endl;
	cout << "l:" << l << endl;
	long long int edge_number = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < G[i].size(); j++) {
			cout << "�}" << edge_number + 1 << "�{��:" << i << " " << G[i][j] << endl;
			edge_number++;
		}
	}
	for (int i = 0; i < n; i++) {
		cout << T[i] << endl;
	}
	���͊m�F�I��*/
	/*�������m�F�J�n
	calculate_Monomials();
	for (int i = 0; i < Monomials.size(); i++) {
		cout << Monomials[i] << endl;
	}
	�������m�F�I��*/

	clock_t start = clock();    //���ԑ���J�n

	clock_t end = clock();     // ���ԑ���I��
	cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";

	return 0;
}