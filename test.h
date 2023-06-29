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
			
			int determinant(int ax, int ay, int az, int bx, int by, int bz, int cx, int cy, int cz){
        		return (ax*by*cz)+(ay*bz*cx)+(az*bx*cy)-(cx*by*az)-(bx*ay*cz)-(ax*cy*bz);
			}
			
			int dot(int ax, int ay, int az, int bx, int by, int bz){
				return ax*bx+ay*by+az*bz;
			}
			
			std::vector<int> cross(int ax, int ay, int az, int bx, int by, int bz){
				TwoDimension twod;
				int x = twod.determinant(ay, az, by, bz);
				int y = twod.determinant(ax, az, bx, bz)*-1;
				int z = twod.determinant(ax, ay, bx, by);
				return {x, y, z};
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
		
			int gcd(int a, int b){
				return b==0?a:gcd(b, a%b);
			}
			
			int lcm(int a, int b){
				return a*b/gcd(a, b);
			}
			
			void equation_multiply(vector<int> &v, int r){
				for(auto &i:v) i*=r;
			}
			
			void row_operation(vector<int> a, vector<int> &b, int r){ 
				for(int i = 0; i < b.size(); i++){
					b[i] += a[i]*r;
				}
			}
				
			vector<double> GaussianElimination(vector<vector<int>> &v){ // Nx(N+1)
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
				vector<double> ans;
				for(int i = 0; i < v.size(); i++){
					double a = v[i][v.size()]/v[i][i];
					ans.emplace_back(a);
				}
				return ans;
			}
			
    };

}

#endif // test
