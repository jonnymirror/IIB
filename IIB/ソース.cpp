#include <iostream>
#include <vector>
#include <algorithm>
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
long long int a=5;//mod�̌v�Z�ɗp����萔
vector<long long int> Monomials;//�n�b�V���l�̌v�Z�ɗp����P����(Monomials[i]��x^i��\��)
vector<long long int> Hash_table;//polynominal_hash_function�ŋ��߂��n�b�V���l(���ꂼ��̒��_�̋ߖT�ɑ΂���)���i�[����\
vector<long long int> Place_of_vertices;//���_i��������^�C�v�p�[�e�B�V������\���z��,Place_of_vertices[i]=-1�̎��͂��̒��_�͂܂��ǂ̃^�C�v�p�[�e�B�V�����ɂ������Ă��Ȃ����Ƃ�\��,Place_of_vertices[i]=j�̎��͒��_i���^�C�v�p�[�e�B�V����j�ɑ����Ă��邱�Ƃ�\��
vector<vector<long long int>> Type_partitions;//�^�C�v�p�[�e�B�V����
long long int nd = -1;//�O���t�̋ߖT���l���̒l���i�[����ϐ�
long long int numbers_of_induced_subgraph_vertices;//�U�������O���t�̒��_��
vector<vector<long long int>> Induced_subgraph;//�U�������O���t
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
void polynominal_hash_fanction() {
	for (long long int i = 0; i < n; i++) {
		long long int sum = 0;
		for (long long int j = 0; j < G[i].size(); j++) {
			sum = (sum + Monomials[G[i][j]]) % mod;
		}
		Hash_table.push_back(sum);
	}
}

//�ߖT���l�������߂�
	void calculate_neighborhood_diversity() {
	long long int count = 0;//�V�����^�C�v�̒��_��type_partition[count]�ɑ�����
	for (long long int v = 0; v < n-1; v++) {
		for (long long int u = v+1; u < n; u++) {
			long long int one = Hash_table[u] - Monomials[v];//h(N(u)\v)
			long long int the_other = Hash_table[v] - Monomials[u];//h(N(v)\u)
			if (one < 0)one += mod;
			if (the_other < 0)the_other += mod;
			one = one % mod;
			the_other = the_other % mod;
			if (one == the_other) {//h(N(u)\v)==h(N(v)\u)�̎�
				if (Place_of_vertices[u] == -1) {//���_u���ǂ̃^�C�v�p�[�e�B�V�����ɂ������Ă��Ȃ���
					Place_of_vertices[u] = count;
					count++;
				}
				Place_of_vertices[v] = Place_of_vertices[u];//v��u�Ɠ����^�C�v�p�[�e�B�V�����ɑ�����
			}
		}
	}
	nd = count;//count�̒l�͋ߖT���l���Ȃ̂�,nd��count�̒l���i�[����
	//�^�C�v�p�[�e�B�V�����̍쐬�J�n
	Type_partitions.resize(count);
	for (long long int i = 0; i < n; i++) {
		Type_partitions[Place_of_vertices[i]].push_back(i);
	}
	//�^�C�v�p�[�e�B�V�����̍쐬�I��
}
//2�̒��_u,v�̋ߖT���r����֐�(�O���t�̗אڃ��X�g���l�̏���������sort���Ă���O��őO�����r���Ă���)
	bool check_neighbor(long long int u, long long int v) {
		long long int count_u = 0;//G[u]�̉��Ԗڂ̒��_����\��
		long long int count_v = 0;//G[v]�̉��Ԗڂ̒��_����\��
		while (1) {
			if (count_u == G[u].size() && count_v == G[v].size())return true;//2�̒��_u,v�̋ߖT�ɑ����钸�_���S�ē�����������,true��Ԃ�
			if (G[u][count_u] == v)count_u++;//u�̋ߖT�ɑ�����v�͔�΂�
			if (G[v][count_v] == u)count_v++;//v�̋ߖT�ɑ�����u�͔�΂�
			if (G[u][count_u] == G[v][count_v]) {
				count_u++;
				count_v++;
			}
			else {
				return false;
			}
		}
}

//�ߖT���l���������Ă��邩�ǂ����m�F����
bool check_neighborhood_diversity() {
	long long int representative = -1;//�����ׂ����^�C�v�p�[�e�B�V�����ɑ������\�_
	long long int others = -1;//�����ׂ����^�C�v�p�[�e�B�V�����ɑ������\�_�ȊO�̒��_
	//���_representative�ƒ��_others�̋ߖT����v���Ă��邩�ǂ����m�F����(�O���t�̗אڃ��X�g���l�̏���������sort���Ă���O��őO�����r���Ă���)
	for (long long int i = 0; i < Type_partitions.size(); i++) {
		if (Type_partitions[i].size() == 1)continue;
		representative = Type_partitions[i][0];//��\�_����܂�
		for (long long int j = 1; j < Type_partitions[i].size(); j++) {
			others = Type_partitions[i][j];//�ߖT���r���钸�_����܂�
			//2�̒��_�̋ߖT�̔�r�����ۂɍs��
			if (!check_neighbor(representative, others)) {//�����قȂ��Ă�����
				return false;
			}
		}
	}
	return true;
}

//�ߖT���l�������߂�,�m�F����,�������ߖT���l�������܂�Ȃ�������x��a��ς��ċ��ߒ���(�ߖT���l���Ɋւ���֐����܂Ƃ߂�����)
void summarize_neighbor_diversity() {
	calculate_Monomials();
	polynominal_hash_fanction();
	calculate_neighborhood_diversity();
	while (!check_neighborhood_diversity()) {//�ߖT���l�����������Ȃ�������ߖT���l�������߂Ȃ���
		//�ϐ���������Ԃɖ߂�(Monomials,Hash_table,Place_of_vertices,Type_partitions,nd��calculate_Monomials�����s����O�ɖ߂�)�J�n
		Monomials.clear();
		Monomials.shrink_to_fit();
		Hash_table.clear();
		Hash_table.shrink_to_fit();
		Place_of_vertices.clear();
		Place_of_vertices.shrink_to_fit();
		for (int i = 0; i < Type_partitions.size(); i++) {
			Type_partitions[i].clear();
			Type_partitions[i].shrink_to_fit();
		}
		nd = -1;
		//�ϐ���������Ԃɖ߂��I��
		//�ߖT���l�������߂Ȃ����J�n
		x++;
		a++;
		calculate_Monomials();
		polynominal_hash_fanction();
		calculate_neighborhood_diversity();
		//�ߖT���l�������߂Ȃ����I��
	}
}

//�U�������O���t���v�Z����֐�
vector<bool> making_induced_subgraph(vector<long long int> Vertex_Subset) {
	vector<bool> Exist;//�U�������O���t���v�Z����̂ɗp����z��,���_i���U�������O���t�̒��_�W���ɑ����Ă�����exist[i]=true,�����Ă��Ȃ����exist[i]=false�ł���.
	Exist.resize(n);
	cout << "Exist�̑傫����" << Exist.size() << endl;
	for (long long int i = 0; i < n; i++) {
		Exist[i] = false;
	}
	for (long long int i = 0; i < Vertex_Subset.size(); i++) {
		Exist[i] = true;
	}
	Induced_subgraph.resize(n);
	for (long long int i = 0; i < n; i++) {
		for (long long int j = 0; j < G[i].size(); j++) {
			//i��G[i][j]���ǂ�����U�������O���t�̒��_�W���ɑ�����Ȃ��
			if (Exist[i] && Exist[G[i][j]]) {
				Induced_subgraph[i].push_back(G[i][j]);
			};//
		}
	}
	return Exist;
}

//�ŏI�I�Ɋ������钸�_�����߂�֐�(bit���g��Ȃ�����)
vector<bool> who_is_influenced_not_bit(vector<bool> Exist) {
	vector<bool> Influenced;//���_i���������Ă�����influenced[i]=true,�����łȂ����influenced[i]=false�ł���
	Influenced.resize(n);
	for (long long int i = 0; i < n; i++) {
		Influenced[i] = false;
	}
	//1��ڂ̊g�U�ߒ��̎����J�n
	for (long long int i = 0; i < n; i++) {
		if (Exist[i] && T[i] == 0)Influenced[i] = true;
	}
	//1��ڂ̊g�U�ߒ��̎����I��
	//t��ڂ̊g�U�ߒ��̎����J�n
	bool changed = false;
	do {
		changed = false;
		for (long long int i = 0; i < n; i++) {
			long long int count = 0;//���_u�̗אړ_�Ŋ������Ă��钸�_�̐��𐔂���
			if (Exist[i] && !Influenced[i]) {
				for (long long int j = 0; j < Induced_subgraph[i].size(); j++) {
					if (Influenced[Induced_subgraph[i][j]])count++;
				}
				if (count >= T[i]) {
					Influenced[i] = true;
					changed = true;
				}
			}
		}
	} while (changed);
	//t��ڂ̊g�U�ߒ��̎����I��
	return Influenced;
}

//Y(X)�����߂�֐�
vector<bool> calculate_YX(vector<bool> Influenced, vector<bool> Exist) {
	vector<bool> YX;//���_i��Y(X)�ɑ�����Ȃ�YX[i]=1,�����Ȃ��Ȃ�YX[i]=0
	YX.resize(n);
	for (long long int i = 0; i < n; i++) {
		YX[i] = false;
	}
	for (long long int i = 0; i < n; i++) {
		if (!Exist[i]) {
			long long int count = 0;
			for (long long int j = 0; j < G[i].size(); j++) {
				if (Influenced[G[i][j]])count++;
			}
			if (count >= T[i])YX[i] = true;
		}
	}
	return YX;
}

//IIB_k�̑O����(G�̃^�C�v�p�[�e�B�V����{V_0,V_1,...,V_nd}�̂��ꂼ���V_i={v_{i,1},...,v_{i,|V_i|}}�̒��_��臒l�̔񌸏���,�Ⴆ��,t(v_{i,j})<=t(v_{i,j+1})),�̂悤�ɕ��ׂ�)
void sort_in_order_of_thresholds() {
	for (long long int i = 0; i < Type_partitions.size(); i++) {
		sort(Type_partitions[i].begin(), Type_partitions[i].end());
	}
}

//�d���g�ݍ��킹��
vector<long long int> A;
bool overlapping_combination(long long int s, long long int t) {
	if (s == nd && t == 0) {
		vector<long long int> Vertex_Subset;
		for (long long int i = 0; i < nd; i++) {
			for (long long int j = 0; j < A[i]; i++) {
				Vertex_Subset.push_back(Type_partitions[i][j]);
			}
		}
		vector<bool> Exist;
		Exist = making_induced_subgraph(Vertex_Subset);
		vector<bool> Influenced;
		Influenced = who_is_influenced_not_bit(Exist);
		vector<bool> YX;
		YX = calculate_YX(Influenced,Exist);
		long long int count = 0;
		for (long long int i = 0; i < n; i++) {
			if (YX[i])count++;
		}
		if (count <= l)return true;
		return false;
	}
	if (s == nd && t != 0) {
		return false;
	}
	for (long long int i = 0; i <= t; i++) {
		A[s] = i;
		overlapping_combination(s + 1, t - i);
	}
}

//IIB_k(G,k,l)//G,k,l�̓O���[�o���ϐ��Őݒ肵�Ă���̂�,�֐��̈����ɏ����Ă��Ȃ�
bool IIB_k() {
	sort_in_order_of_thresholds();
	for (long long int f = 1; f < k + 1; f++) {
		if (overlapping_combination(0, f))return true;
	}
	return false;
}

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
	mod =a* n * n;//mod�͒萔�~n^2

	for(int i = 0; i < n; i++){
		Place_of_vertices.push_back(-1);//
	}
	//���͏I��

	//���͐����J�n
	for (int i = 0; i < n; i++) {
		sort(G[i].begin(), G[i].end());//�O���t�̗אڃ��X�g��l�����������ɕ��ׂ�(�ߖT���l���������Ă��邩�ǂ����̊m�F�����鎞�̂��߂�)
	}
	//���͐����I��
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
	for (int i = 0; i < Monomials.size(); i++) {
		cout << Monomials[i] << endl;
	}
	�������m�F�I��*/

	clock_t start = clock();    //���ԑ���J�n
	
	if (IIB_k) {
		cout << "YES" << endl;
	}
	else {
		cout << "No" << endl;
	}

	clock_t end = clock();     // ���ԑ���I��
	cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";

	return 0;
}