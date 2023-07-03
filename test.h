#include<iostream>
#include<vector>
#include<math.h>
#ifndef test
#define test
#define int long long

using namespace std;

namespace lib{

    class TwoDimension{

        public:
        	
        	int determinant(int ax, int ay, int bx, int by){
        		return ax*by-ay*bx;
			}
            
            double dis_from_point_to_point(int ax, int ay, int bx, int by){
                double dis = sqrt(1.0*((ax-bx)*(ax-bx) + (ay-by)*(ay-by))); 
                return dis;
            }

            double dis_from_point_to_line(int a, int b, int c, int x, int y){   //ax+by+c=0, p(x, y)
                double mom = sqrt(1.0*(a*a+b*b));
                double son = abs(a*x+b*y+c);
                return son/mom;
            }

    };

    class ThreeDimension{

        public:
			
			int dot(int ax, int ay, int az, int bx, int by, int bz){
				return ax*bx+ay*by+az*bz;
			}
			
			vector<int> cross(int ax, int ay, int az, int bx, int by, int bz){
				TwoDimension twod;
				int x = twod.determinant(ay, az, by, bz);
				int y = twod.determinant(ax, az, bx, bz)*-1;
				int z = twod.determinant(ax, ay, bx, by);
				return {x, y, z};
			}
			
			int determinant(int ax, int ay, int az, int bx, int by, int bz, int cx, int cy, int cz){
        		vector<int> v = cross(ax, ay, az, bx, by, bz);
        		return dot(v[0], v[1], v[2], cx, cy, cz);
			}
			
            double dis_from_point_to_point(int ax, int ay, int az, int bx, int by, int bz){
                double dis = sqrt(1.0*((ax-bx)*(ax-bx) + (ay-by)*(ay-by) + (az-bz)*(az-bz))); 
                return dis;
            }

            double dis_from_point_to_plane(int a, int b, int c, int d, int x, int y, int z){   //ax+by+cz+d=0, p(x, y, z)
                double mom = sqrt(1.0*(a*a+b*b+c*c));
                double son = abs(a*x+b*y+c*z+d);	
                return son/mom;
            }
    
    };


    class Matrix{
		
		public:
			
			vector<vector<int>> unitmatrix(int n){
				vector<vector<int>> v(n);
				for(int i = 0; i < n; i++){
					for(int j = 0; j < n; j++){
						if(i==j) v[i].push_back(1);
						else v[i].push_back(0);
					}
				}
				return v;
			}
		
			int gcd(int a, int b){
				return b==0?a:gcd(b, a%b);
			}
			
			void make_equation_simple(vector<int> &v){
				int g = gcd(v[0], v[1]);
				for(int i = 2; i < v.size(); i++){
					g = gcd(g, v[i]);
				}
				for(auto &i:v) i/=g;
			}
			
			int lcm(int a, int b){
				return a*b/gcd(a, b);
			}
			
			void equation_multiply(vector<int> &v, int k){
				for(auto &i:v) i*=k;
			}
			
			void row_operation(vector<int> a, vector<int> &b, int k){ 
				for(int i = 0; i < b.size(); i++){
					b[i] += a[i]*k;
				}
				make_equation_simple(b);
			}
				
			vector<double> GaussianElimination(vector<vector<int>> &v){ // Nx(N+1)
				vector<double> ans;
				for(int i = 0; i < v.size()-1; i++){
					for(int j = i+1; j < v.size(); j++){
						if(v[i][i]<0) equation_multiply(v[i], -1);
						if(v[j][i]<0) equation_multiply(v[j], -1);
						
						if(v[j][i]%v[i][i]==0){
							row_operation(v[i], v[j], v[j][i]/v[i][i]*-1);
							continue;
						}
						int t = lcm(v[i][i], v[j][i]);
						equation_multiply(v[j], t/v[j][i]);
						row_operation(v[i], v[j], v[j][i]/v[i][i]*-1);
					}
				}
				
				
				if(v[v.size()-1][v.size()]==0&&v[v.size()-1][v.size()-1]==0){
					cout << "infinite solutions\n";
					return ans;
				}
				
				if(v[v.size()-1][v.size()]!=0&&v[v.size()-1][v.size()-1]==0){
					cout << "no solution\n";
					return ans;
				}
				
				for(int i = v.size()-1; i >= 1; i--){
					for(int j = i-1; j >=0; j--){
						
						if(v[i][i]<0)
						{
							equation_multiply(v[i], -1);
						 } 
						if(v[j][i]<0){
							equation_multiply(v[j], -1);
						} 
						
						
						if(v[j][i]%v[i][i]==0){
							row_operation(v[i], v[j], v[j][i]/v[i][i]*-1);
							continue;
						}
						
						int t = lcm(v[i][i], v[j][i]);
						equation_multiply(v[j], t/v[j][i]);
						row_operation(v[i], v[j], v[j][i]/v[i][i]*-1);
					}
				}
				
				
				
				for(int i = 0; i < v.size(); i++){
					double a = v[i][v.size()]/v[i][i];
					ans.emplace_back(a);
				}
				return ans;
				
			}
			
			vector<vector<int>> matrix_multiply(vector<vector<int>> a, vector<vector<int>> b){
				vector<vector<int>> c(a.size());
				for(int i = 0; i < a.size(); i++){
					for(int j = 0; j < a.size(); j++){
						int tmp = 0;
						for(int k = 0; k < b.size(); k++){
							tmp += a[i][k]*b[k][j];
						}
						c[i].emplace_back(tmp);
					}
				}
				return c;
			}
			
			vector<vector<int>> matrix_fast_pow(vector<vector<int>> a, int n){
				if(n==0) return unitmatrix(a.size());
				if(n==1) return a;
				vector<vector<int>> tmp = matrix_multiply(a, a);
				if(n%2==0) return matrix_fast_pow(tmp, n/2);
				return matrix_multiply(matrix_fast_pow(tmp, (n-1)/2), a);
			}
    };

}

#endif // test
