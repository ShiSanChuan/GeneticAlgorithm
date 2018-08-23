#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include "opencv2/opencv.hpp"
#define pi 3.1415926

class GA
{
private:
	enum rank_menthod{UP=0,DOWN=1};
	int chrom_num;
	int gene_num;
	float p_recombin;
	float p_mut;
	float search_min;
	float search_max;
	std::function<float (float&)> fun;
public:
	GA(std::function<float (float&)> _fun,const int _chrom_num=40,const int _gene_num=20,
		const float _p_recombin=0.1,const float _p_mut=0.2,
		const float min=0,const float max=10);
	~GA(){};
	cv::Mat crtbp(int Nind=0,int Lind=0);
	std::pair<float, float> ranking(std::vector<float> &objV,rank_menthod method=UP);
	void select(cv::Mat &Popula,std::vector<float> &rank);
	void recombin(cv::Mat &Popula,float opt=0);
	void mut(cv::Mat &Popula,float opt=0);
	std::vector<float> bs2rv(cv::Mat &Popula,float min=0,float max=0);
};