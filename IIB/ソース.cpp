#include <iostream>
#include <vector>
#include <intrin.h>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>  

using namespace std;

//定数開始
vector<vector<long long int>> G;//入力されるグラフ
long long int n;//頂点数
long long int m;//枝数
long long int k;//パラメータ(最終的に感染する頂点の数の最大値)
long long int l;//パラメータ(免疫を持たせる頂点の数の最大値)
vector<long long int> T;//関数
long long int x = 5;//ハッシュ値の計算に用いるhash-multiplier
long long int mod;//ハッシュ値の計算に用いるmod(頂点数によって変化する)
vector<long long int> Monomials;//ハッシュ値の計算に用いる単項式(Monomials[i]はx^iを表す)
//定数終了
//最終的に感染する頂点を求める関数
long long int who_is_influenced(long long int bit) {
	long long int influenced = 0;//感染している頂点をbit1で表す
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

//Monomialsを求める
void calculate_Monomials() {
	long long int count = 1;
	for (long long int i = 0; i < n; i++) {
		Monomials.push_back(count);
		count = (count*x)%mod;
	}
}
//polynominal hash function(それぞれの頂点の近傍に対するハッシュ値を求める)


//メイン関数

int main() {
	//入力開始
	ifstream ifs("paper_sample.txt");

	if (!ifs) {
		std::cout << "Error!";
		return 1;
	}

	string s;
	long long int count = 0;
	long long int x = 0;
	long long int y = 0;

	while (getline(ifs, s, ' ')) {     // スペース（' '）で区切って，格納
		cout << "count:" << count << endl;
		if (count == 0) {//頂点数
			n = stoll(s);
			G.resize(n);//グラフの大きさ確保
			count++;
		}
		else if (count == 1) {//枝数
			m = stoll(s);
			count++;
		}
		else if (count == 2) {//パラメータk
			k = stoll(s);
			count++;
		}
		else if (count == 3) {//パラメータl
			l = stoll(s);
			count++;
		}
		else if (count > 3 && count < 4 + 2 * m && count % 2 == 0) {//グラフの枝の端点
			x = stoll(s);
			x--;
			count++;
		}
		else if (count > 3 && count < 4 + 2 * m && count % 2 == 1) {//グラフの枝のもう一つの端点
			y = stoll(s);
			y--;
			count++;
			G[x].push_back(y);
			G[y].push_back(x);
		}
		else {//頂点の閾値
			x = stoll(s);
			T.push_back(x);
			count++;
		}
	}
	mod = n * n;//modは定数×n^2
	//入力終了

	/*入力確認開始
	cout << "頂点数:" << n << endl;
	cout << "枝数:" << m << endl;
	cout << "k:" << k << endl;
	cout << "l:" << l << endl;
	long long int edge_number = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < G[i].size(); j++) {
			cout << "枝" << edge_number + 1 << "本目:" << i << " " << G[i][j] << endl;
			edge_number++;
		}
	}
	for (int i = 0; i < n; i++) {
		cout << T[i] << endl;
	}
	入力確認終了*/
	/*多項式確認開始
	calculate_Monomials();
	for (int i = 0; i < Monomials.size(); i++) {
		cout << Monomials[i] << endl;
	}
	多項式確認終了*/

	clock_t start = clock();    //時間測定開始

	clock_t end = clock();     // 時間測定終了
	cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";

	return 0;
}