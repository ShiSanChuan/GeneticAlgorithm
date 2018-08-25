#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include "opencv2/opencv.hpp"
#define pi 3.1415926

class GA
{
protected:
	int chrom_num;
	int gene_num;
	float p_recombin;
	float p_mut;
	float search_min;
	float search_max;
	float (*fun)(std::vector<float> argv);
	int para_num;
	cv::Mat ost;//fun ,x,y,z....
public:
	GA(const int _chrom_num=40,const int _gene_num=20,
		const float _p_recombin=0.3,const float _p_mut=0.2,
		const float min=0,const float max=1,const int _para_num=1);
	void solve(float (*_fun)(std::vector<float> argv),const int &_para_num=0);
	~GA(){};
	cv::Mat crtbp(const int &Nind=0,const int &Lind=0);
	std::pair<std::vector<float> , float> ranking(void);
	GA& select(cv::Mat &Popula);
	GA& recombin(cv::Mat &Popula,const float &opt=0);
	GA& mut(cv::Mat &Popula,float opt=0);
	void bs2rv(cv::Mat &Popula,float min=0,float max=0);
};

class GA_BP:public GA
{
private:
	cv::Mat input;
	cv::Mat output;
	int implication_num;
public:
	GA_BP(const int _chrom_num=40,const int _gene_num=20,
		const float _p_recombin=0.3,const float _p_mut=0.2,
		const float min=0,const float max=1,const int _para_num=1);
	~GA_BP(){};
	void BPsolve(cv::Mat &_input,cv::Mat &_output);
	std::pair<std::vector<float>, float> ranking(void);
};