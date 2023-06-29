#include <iostream>
#include "test.h"
#define int long long
using namespace std;
using namespace lib;


signed main() {
  TwoDimension twod;
  ThreeDimension threed;
  Matrix mtx;
  int n; cin >> n;
  vector<vector<int>> v(n);
  for(int i = 0; i<n; i++){
	for(int j = 0; j < n+1; j++){
  		int a; cin >> a;
  		v[i].emplace_back(a);
		}
	}
	vector<double> ans = mtx.GaussianElimination(v);
	for(auto i:ans){
		cout << i << '\n';
	}
  return 0;
}
