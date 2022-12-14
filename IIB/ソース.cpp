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
//最終的に感染する頂点を求める関数(使っていない)
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
	cout << "calculate_Monomials開始" << endl;
	long long int count = 1;
	for (long long int i = 0; i < n; i++) {
		cout <<"i="<<i << endl;
		Monomials.push_back(count);
		count = (count*x)%mod;
	}
	cout << "calculate_Monomials終了" << endl;
}

//polynominal hash function(それぞれの頂点の近傍に対するハッシュ値を求める)
void polynominal_hash_fanction() {
	cout << "polynominal_hash_function開始" << endl;
	for (long long int i = 0; i < n; i++) {
		cout <<"i="<<i << endl;
		long long int sum = 0;
		for (long long int j = 0; j < G[i].size(); j++) {
			cout << "j=" << j << endl;
			sum = (sum + Monomials[G[i][j]]) % mod;
		}
		Hash_table.push_back(sum);
	}
	cout << "polynominal_hash_function終了" << endl;
}

//近傍多様性を求める
void calculate_neighborhood_diversity() {
	cout << "calculate_neighborhood_diversity開始" << endl;
	long long int count = 0;//新しいタイプの頂点はtype_partition[count]に属する
	Place_of_vertices[0] = 0;
	for (long long int v = 0; v < n; v++) {
		cout << "v=" << v << endl;
		for (long long int u = 0; u < n; u++) {
			cout << "u=" << u << endl;
			long long int one = Hash_table[u] - Monomials[v];//h(N(u)\v)
			long long int the_other = Hash_table[v] - Monomials[u];//h(N(v)\u)
			if (one < 0)one += mod;
			if (the_other < 0)the_other += mod;
			one = one % mod;
			the_other = the_other % mod;
			cout << "one=" << one << endl;
			cout << "the_other=" << the_other << endl;
			if (one == the_other) {//h(N(u)\v)==h(N(v)\u)の時
				cout << "Place_of_vertices[" << v << "]=" << Place_of_vertices[v]<<endl;
				cout << "Place_of_vertices[" << u << "]=" << Place_of_vertices[u]<<endl;
				if(Place_of_vertices[u]!=-1){
					Place_of_vertices[v] = Place_of_vertices[u];//vはuと同じタイプパーティションに属する
				}
				cout << "changed!" << endl;
				cout << "Place_of_vertices[" << v << "]=" << Place_of_vertices[v] << endl;
				cout << "Place_of_vertices[" << u << "]=" << Place_of_vertices[u] << endl;
			}
		}
		if (Place_of_vertices[v] == -1) {//頂点vがどのタイプパーティションにも属していない時
			count++;
			cout << "count=" << count << endl;
			Place_of_vertices[v] = count;
		}
	}
	nd = count+1;//countの値は近傍多様性なので,ndにcountの値を格納する
	//Place_of_verticesの確認
	for (long long int i = 0; i < n; i++) {
		cout << "Place_of_vertices[" << i << "]="<<Place_of_vertices[i] << endl;
	}
	//タイプパーティションの作成開始
	Type_partitions.resize(nd);
	for (long long int i = 0; i < n; i++) {
		cout << "i=" << i << endl;
		Type_partitions[Place_of_vertices[i]].push_back(i);
	}
	//タイプパーティションの作成終了
	cout << "calculate_neighborhood_diversity終了" << endl;
}
//2つの頂点u,vの近傍を比較する関数(グラフの隣接リストが値の小さい順にsortしている前提で前から比較していく)
bool check_neighbor(long long int u, long long int v) {
	cout << "check_neighbor開始" << endl;
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
	cout << "check_neighbor終了" << endl;
}

//近傍多様性が合っているかどうか確認する
bool check_neighborhood_diversity() {
	cout << "check_neighborhood_diversity開始" << endl;
	long long int representative = -1;//今調べたいタイプパーティションに属する代表点
	long long int others = -1;//今調べたいタイプパーティションに属する代表点以外の頂点
	//頂点representativeと頂点othersの近傍が一致しているかどうか確認する(グラフの隣接リストが値の小さい順にsortしている前提で前から比較していく)
	for (long long int i = 0; i < Type_partitions.size(); i++) {
		cout << "i=" << i << endl;
		if (Type_partitions[i].size() == 1)continue;
		representative = Type_partitions[i][0];//代表点が定まる
		for (long long int j = 1; j < Type_partitions[i].size(); j++) {
			cout << "j=" << j <<endl;
			others = Type_partitions[i][j];//近傍を比較する頂点が定まる
			//2つの頂点の近傍の比較を実際に行う
			if (!check_neighbor(representative, others)) {//もし異なっていたら
				cout << "異なっているぞ！" <<endl;
				cout << "check_neighborhood_diversity終了" << endl;
				return false;
			}
		}
	}
	cout << "check_neighborhood_diversity終了" << endl;
	return true;
}

//近傍多様性を求めて,確認する,正しい近傍多様性が求まらなかったらxやaを変えて求め直す(近傍多様性に関する関数をまとめたもの)
void summarize_neighbor_diversity() {
	cout << "summarize_neighbor_diversity開始" << endl;
	calculate_Monomials();
	polynominal_hash_fanction();
	calculate_neighborhood_diversity();
	while (!check_neighborhood_diversity()) {//近傍多様性が正しくなかったら近傍多様性を求めなおす
		cout << "近傍多様性が正しくないのだ！" << endl;
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
	cout << "summarize_neighbor_diversity終了" << endl;
}

//誘導部分グラフを計算する関数
vector<bool> making_induced_subgraph(vector<long long int> Vertex_Subset) {
	vector<bool> Exist;//誘導部分グラフを計算するのに用いる配列,頂点iが誘導部分グラフの頂点集合に属していたらexist[i]=true,属していなければexist[i]=falseである.
	Exist.resize(n);
	for (long long int i = 0; i < n; i++) {
		Exist[i] = false;
	}
	for (long long int i = 0; i < Vertex_Subset.size(); i++) {
		Exist[Vertex_Subset[i]] = true;
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
	return Exist;
}

//最終的に感染する頂点を求める関数(bitを使わない実装)
vector<bool> who_is_influenced_not_bit(vector<bool> Exist) {
	vector<bool> Influenced;//頂点iが感染していたらinfluenced[i]=true,そうでなければinfluenced[i]=falseである
	Influenced.resize(n);
	for (long long int i = 0; i < n; i++) {
		Influenced[i] = false;
	}
	//1回目の拡散過程の実装開始
	for (long long int i = 0; i < n; i++) {
		if (Exist[i] && T[i] == 0)Influenced[i] = true;
	}
	//1回目の拡散過程の実装終了
	//t回目の拡散過程の実装開始
	bool changed = false;
	do {
		changed = false;
		for (long long int i = 0; i < n; i++) {
			long long int count = 0;//頂点uの隣接点で感染している頂点の数を数える
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
	//t回目の拡散過程の実装終了
	return Influenced;
}

//Y(X)を求める関数
vector<bool> calculate_YX(vector<bool> Influenced, vector<bool> Exist) {
	vector<bool> YX;//頂点iがY(X)に属するならYX[i]=1,属さないならYX[i]=0
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

//IIB_kの前処理(Gのタイプパーティション{V_0,V_1,...,V_nd}のそれぞれのV_i={v_{i,1},...,v_{i,|V_i|}}の頂点を閾値の非減少順,例えば,t(v_{i,j})<=t(v_{i,j+1})),のように並べる)
void sort_in_order_of_thresholds() {
	cout << "sort_in_order_of_thresholds開始" << endl;
	for (long long int i = 0; i < Type_partitions.size(); i++) {
		sort(Type_partitions[i].begin(), Type_partitions[i].end());
	}
	cout << "sort_in_order_of_thresholds終了" <<endl;
}

//重複組み合わせ列挙
vector<long long int> A;
bool overlapping_combination(long long int s, long long int t) {
	cout << "overlapping_combination 開始　s=" << s << " t=" << t << endl;
	if (s == nd && t == 0) {
		//デバック用開始
		for (long long int i = 0; i < nd; i++) {
			cout << A[i] << " ";
		}
		cout << endl;
		//デバック用終了
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
		//デバック開始
		cout << "Vertex_Subset" << endl;
		for (long long int i = 0; i < Vertex_Subset.size(); i++) {
			cout << Vertex_Subset[i]+1 << " ";
		}
		cout << endl;
		cout << "Exist" << endl;
		for (long long int i = 0; i < n; i++) {
			cout << Exist[i] << " ";
		}
		cout << endl;
		cout << "Influenced" << endl;
		for (long long int i = 0; i < n; i++) {
			cout << Influenced[i] << " ";
		}
		cout << endl;
		cout << "YX" << endl;
		for (long long int i = 0; i < n; i++) {
			cout << YX[i] << " ";
		}
		cout << endl;
		cout << "count=" << count << endl;
		cout << "l=" << l << endl;
		cout << "k=" << k << endl;
		//デバック終了
		if (count <= l) {
			cout << "trueです" << endl;
			return true;
		}
		cout << "falseです.aaaaaaaaaa" << endl;
		return false;
	}
	if (s == nd && t != 0) {
		cout << "falseです.bbbbbbbbbb" << endl;
		return false;
	}
	for (long long int i = 0; i <= t; i++) {
		cout << "i=" << i << endl;
		A[s] = i;
		return overlapping_combination(s + 1, t - i);
	}
}

//IIB_k(G,k,l)//G,k,lはグローバル変数で設定しているので,関数の引数に書いていない
bool IIB_k() {
	cout << "IIB_k開始" << endl;
	sort_in_order_of_thresholds();
	for (long long int f = 1; f < k + 1; f++) {
		cout << "f=" << f << endl;
		if (overlapping_combination(0, f)) {
			cout << "trueでIIB_k終了" << endl;
			return true;
		}
	}
	cout << "falseでIIB_k終了" << endl;
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
	
	summarize_neighbor_diversity();

	A.resize(n);

	if (IIB_k()) {
		cout << "YES" << endl;
	}
	else {
		cout << "No" << endl;
	}

	clock_t end = clock();     // 時間測定終了
	cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";

	return 0;
}