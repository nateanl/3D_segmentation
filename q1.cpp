#include<iostream>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include<vector>
#include <Eigen/Dense>
#include <fstream>
using namespace std;
using namespace Eigen;

typedef unordered_map<pair<int,int>,Vector3d,boost::hash<pair<int, int>>> points;
typedef unordered_map<pair<int,int>,bool,boost::hash<pair<int, int>>> empty;
pair<int,int> get_col_row(string file_dir){
	ifstream f;
	f.open(file_dir);
	int col=0, row =0;
	f>>col;
	f>>row;
	f.close();
	return make_pair(col,row);
}
points load_data(string file_dir){
	points data;
	//open the file
	ifstream f;
	f.open(file_dir);
	string line;
	int col=0, row =0;
	f>>col;
	f>>row;
	cout<<col<<" "<<row<<endl;
	for(int i =0; i<9; i++){
		getline(f,line);
	}
	double x=0.0,y=0.0,z=0.0,intensity=0.0;
	for(int i = 0; i<col; i++){
		for(int j = 0; j<row; j++){
			f>>x>>y>>z>>intensity;
			data[make_pair(i,j)] = Vector3d(x,y,z);
		}
	}
	f.close();
	return data;
}
Vector3d get_normal_vector(list<Vector3d> &neighbors){
	Vector3d center(0.0,0.0,0.0);
	Matrix3d cov;
	for(auto p : neighbors){
		center+=p;
	}
	center = center/neighbors.size();
	for(auto p : neighbors){
		Vector3d pprime = p-center;
		cov += pprime * pprime.transpose();
	}
	EigenSolver<MatrixXd> es(cov);
	Vector3d eigenvalue = es.eigenvalues().real();
	int min_index = 0;
	for(int i = 1; i<eigenvalue.size(); i++){
		if(eigenvalue(i)<eigenvalue(min_index)){
			min_index = i;
		}
	}
	Matrix3d eigenmatrix = es.eigenvectors().real();
	return eigenmatrix.col(min_index);
}
points compute_normal(points data, empty is_empty, int col, int row){
	points normal;
	list<Vector3d> neighbors;
	for(auto i : data){
		int c = i.first.first;
		int r = i.first.second;
		if(!is_empty[i.first] && c-1>=0 && c+1<=col && r-1>=0 && r+1<row){
			for(int j = -1; j<=1; j++){
				for(int k = -1; k<=1; k++){
					if(!is_empty[make_pair(c+j,r+k)]){
						neighbors.push_back(data[make_pair(c+j,r+k)]);
					}
				}
			}
		}
		if(neighbors.size()>1){
			normal[i.first]=get_normal_vector(neighbors);
		}
	}
	return normal;
}

empty check_empty(points data){
	empty is_empty;
	long long count = 0;
	for(auto point : data){
		if(abs(point.second[0])<0.1 && abs(point.second[1])<0.1 && abs(point.second[2])<0.1){
			is_empty[point.first] = true;
		}
		else{
			is_empty[point.first] = false;
			count++;
		}
	}
	cout<<count<<endl;
	return is_empty;
}
int main(){
	string file_dir = "./DATA/big_example.ptx";
	pair<int,int> c_r = get_col_row(file_dir);
	cout<<c_r.first<<" "<<c_r.second<<endl;
	points data = load_data(file_dir);
	empty is_empty = check_empty(data);
	points normal = compute_normal(data,is_empty,c_r.first,c_r.second);
	return 0;
}