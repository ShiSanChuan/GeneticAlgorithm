#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include "opencv2/opencv.hpp"
#define pi 3.1415926
//GA_base
class GA
{
protected:
	int chrom_num;
	int gene_num;
	float p_recombin;
	float p_mut;
	float search_min;
	float search_max;
	int para_num;
	cv::Mat ost;//fun ,x,y,z..
private:
	float (*fun)(std::vector<float> argv);
public:
	GA(const int _chrom_num=40,const int _gene_num=20,
		const float _p_recombin=0.3,const float _p_mut=0.2,
		const float min=0,const float max=1,const int _para_num=1);
	void solve(float (*_fun)(std::vector<float> argv),const int &_para_num=0);
	~GA(){};
	cv::Mat crtbp(const int &Nind=0,const int &Lind=0,const int&encodemin=0,const int &encodemax=2);
	std::pair<std::vector<float> , float> ranking(void);
	GA& select(cv::Mat &Popula,int _method=0);
	GA& recombin(cv::Mat &Popula,const float &opt=0);
	GA& mut(cv::Mat &Popula,float opt=0);
	void bs2rv(cv::Mat &Popula,float min=0,float max=0);
};
//BP
class GA_BP:public GA
{
private:
	cv::Mat input;
	cv::Mat output;
	int implication_num;
public:
	GA_BP(const int _chrom_num=40,const int _gene_num=20,
		const float _p_recombin=0.3,const float _p_mut=0.2,
		const float min=0,const float max=1,const int _para_num=1):
	GA(_chrom_num,_gene_num,_p_recombin,_p_mut,min,max,_para_num){}

	~GA_BP(){};
	void BPsolve(cv::Mat &_input,cv::Mat &_output);
	std::pair<std::vector<float>, float> ranking(void);
};
//GA_TSP
class GA_TSP:public GA
{
private:
	cv::Mat address;
	double distance(int indexi,int indexj);
public:
	GA_TSP(const int _chrom_num=40,const int _gene_num=20,
		const float _p_recombin=0.3,const float _p_mut=0.2,
		const float min=0,const float max=1,const int _para_num=1):
	GA(_chrom_num,_gene_num,_p_recombin,_p_mut,min,max,_para_num){}
	
	~GA_TSP(){};
	void TSPsolve(cv::Mat &_address);
	std::pair<std::vector<float>, float> ranking(cv::Mat &_Poulate);
	cv::Mat crtbp(int encodemax=0);
	GA_TSP& recombin(cv::Mat &Popula,const float &opt=0);
	GA_TSP& mut(cv::Mat &Popula,float opt=0);
	GA_TSP& select(cv::Mat &Popula,int _method=0){GA::select(Popula,_method);return *this;}
};


