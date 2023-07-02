#include <iostream>
#include "test.h"
#define int long long
using namespace std;
using namespace lib;


signed main(){
	ThreeDimension threed;
	TwoDimension twod;
	Matrix mtx; 
	int n; cin >> n;
	vector<vector<int>> v = {{1, 1}, {1, 0}};
	v = mtx.matrix_fast_pow(v, n-1);
	for(int i = 0; i< 2; i++){
		for(int j = 0; j < 2; j++){
			cout << v[i][j] << " ";
		}
		cout << '\n';
	}
}
