#include <iostream>
#include <vector>
#include <algorithm>
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
long long int a=5;//modの計算に用いる定数
vector<long long int> Monomials;//ハッシュ値の計算に用いる単項式(Monomials[i]はx^iを表す)
vector<long long int> Hash_table;//polynominal_hash_functionで求めたハッシュ値(それぞれの頂点の近傍に対する)を格納する表
vector<long long int> Place_of_vertices;//頂点iが属するタイプパーティションを表す配列,Place_of_vertices[i]=-1の時はその頂点はまだどのタイプパーティションにも属していないことを表す,Place_of_vertices[i]=jの時は頂点iがタイプパーティションjに属していることを表す
vector<vector<long long int>> Type_partitions;//タイプパーティション
long long int nd = -1;//グラフの近傍多様性の値を格納する変数
long long int numbers_of_induced_subgraph_vertices;//誘導部分グラフの頂点数
vector<vector<long long int>> Induced_subgraph;//誘導部分グラフ
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
void polynominal_hash_fanction() {
	for (long long int i = 0; i < n; i++) {
		long long int sum = 0;
		for (long long int j = 0; j < G[i].size(); j++) {
			sum = (sum + Monomials[G[i][j]]) % mod;
		}
		Hash_table.push_back(sum);
	}
}

//近傍多様性を求める
	void calculate_neighborhood_diversity() {
	long long int count = 0;//新しいタイプの頂点はtype_partition[count]に属する
	for (long long int v = 0; v < n-1; v++) {
		for (long long int u = v+1; u < n; u++) {
			long long int one = Hash_table[u] - Monomials[v];//h(N(u)\v)
			long long int the_other = Hash_table[v] - Monomials[u];//h(N(v)\u)
			if (one < 0)one += mod;
			if (the_other < 0)the_other += mod;
			one = one % mod;
			the_other = the_other % mod;
			if (one == the_other) {//h(N(u)\v)==h(N(v)\u)の時
				if (Place_of_vertices[u] == -1) {//頂点uがどのタイプパーティションにも属していない時
					Place_of_vertices[u] = count;
					count++;
				}
				Place_of_vertices[v] = Place_of_vertices[u];//vはuと同じタイプパーティションに属する
			}
		}
	}
	nd = count;//countの値は近傍多様性なので,ndにcountの値を格納する
	//タイプパーティションの作成開始
	Type_partitions.resize(count);
	for (long long int i = 0; i < n; i++) {
		Type_partitions[Place_of_vertices[i]].push_back(i);
	}
	//タイプパーティションの作成終了
}
//2つの頂点u,vの近傍を比較する関数(グラフの隣接リストが値の小さい順にsortしている前提で前から比較していく)
	bool check_neighbor(long long int u, long long int v) {
		long long int count_u = 0;//G[u]の何番目の頂点かを表す
		long long int count_v = 0;//G[v]の何番目の頂点かを表す
		while (1) {
			if (count_u == G[u].size() && count_v == G[v].size())return true;//2つの頂点u,vの近傍に属する頂点が全て同じだったら,trueを返す
			if (G[u][count_u] == v)count_u++;//uの近傍に属するvは飛ばす
			if (G[v][count_v] == u)count_v++;//vの近傍に属するuは飛ばす
			if (G[u][count_u] == G[v][count_v]) {
				count_u++;
				count_v++;
			}
			else {
				return false;
			}
		}
}

//近傍多様性が合っているかどうか確認する
bool check_neighborhood_diversity() {
	long long int representative = -1;//今調べたいタイプパーティションに属する代表点
	long long int others = -1;//今調べたいタイプパーティションに属する代表点以外の頂点
	//頂点representativeと頂点othersの近傍が一致しているかどうか確認する(グラフの隣接リストが値の小さい順にsortしている前提で前から比較していく)
	for (long long int i = 0; i < Type_partitions.size(); i++) {
		if (Type_partitions[i].size() == 1)continue;
		representative = Type_partitions[i][0];//代表点が定まる
		for (long long int j = 1; j < Type_partitions[i].size(); j++) {
			others = Type_partitions[i][j];//近傍を比較する頂点が定まる
			//2つの頂点の近傍の比較を実際に行う
			if (!check_neighbor(representative, others)) {//もし異なっていたら
				return false;
			}
		}
	}
	return true;
}

//近傍多様性を求めて,確認する,正しい近傍多様性が求まらなかったらxやaを変えて求め直す(近傍多様性に関する関数をまとめたもの)
void summarize_neighbor_diversity() {
	calculate_Monomials();
	polynominal_hash_fanction();
	calculate_neighborhood_diversity();
	while (!check_neighborhood_diversity()) {//近傍多様性が正しくなかったら近傍多様性を求めなおす
		//変数を初期状態に戻す(Monomials,Hash_table,Place_of_vertices,Type_partitions,ndをcalculate_Monomialsを実行する前に戻す)開始
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
		//変数を初期状態に戻す終了
		//近傍多様性を求めなおす開始
		x++;
		a++;
		calculate_Monomials();
		polynominal_hash_fanction();
		calculate_neighborhood_diversity();
		//近傍多様性を求めなおす終了
	}
}

//誘導部分グラフを計算する関数
void making_induced_subgraph(long long int induced_subgraph_n) {
	vector<bool> Exist;//誘導部分グラフを計算するのに用いる配列,頂点iが誘導部分グラフの頂点集合に属していたらexist[i]=true,属していなければexist[i]=falseである.
	Exist.resize(n);
	for (long long int i = 0; i < n; i++) {
		Exist[i] = false;
	}
	for (long long int i = 0; i < induced_subgraph_n; i++) {
		long long int a;
		cin >> a;
		a--;
		Exist[a] = true;
	}
	Induced_subgraph.resize(n);
	for (long long int i = 0; i < n; i++) {
		for (long long int j = 0; j < G[i].size(); j++) {
			//iとG[i][j]がどちらも誘導部分グラフの頂点集合に属するならば
			if (Exist[i] && Exist[G[i][j]]) {
				Induced_subgraph[i].push_back(G[i][j]);
			};//
		}
	}
}

//IIB_kの前処理(Gのタイプパーティション{V_0,V_1,...,V_nd}のそれぞれのV_i={v_{i,1},...,v_{i,|V_i|}}の頂点を閾値の非減少順,例えば,t(v_{i,j})<=t(v_{i,j+1})),のように並べる)
void sort_in_order_of_thresholds() {

}

//IIB_k(G,k,l)//G,k,lはグローバル変数で設定しているので,関数の引数に書いていない
bool IIB_k() {
	for (long long int f = 1; f < k + 1; f++) {
		//多重組み合わせの列挙方法が分からん
	}
	return false;
}

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
	mod =a* n * n;//modは定数×n^2

	for(int i = 0; i < n; i++){
		Place_of_vertices.push_back(-1);//
	}
	//入力終了

	//入力整理開始
	for (int i = 0; i < n; i++) {
		sort(G[i].begin(), G[i].end());//グラフの隣接リストを値が小さい順に並べる(近傍多様性が合っているかどうかの確認をする時のために)
	}
	//入力整理終了
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
	for (int i = 0; i < Monomials.size(); i++) {
		cout << Monomials[i] << endl;
	}
	多項式確認終了*/

	clock_t start = clock();    //時間測定開始

	cin >> numbers_of_induced_subgraph_vertices;
	making_induced_subgraph(numbers_of_induced_subgraph_vertices);//部分誘導グラフの作成

	clock_t end = clock();     // 時間測定終了
	cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";

	return 0;
}