#include "GA.h"
/**
 * 初始化 init parament
 */
PSO::PSO(const int _chrom_num,const int _para_num,
		const float min,const float max,
		const float _c1,const float _c2,const float _wmax,
		const float _wmin){
	chrom_num=_chrom_num;
	para_num=_para_num;
	c1=_c1;
	c2=_c2;
	wmax=_wmax;
	wmin=_wmin;
	para_num=_para_num;
	bmin=min;
	bmax=max;
}
/**
 * 输入目标参数
 */
void PSO::solve(float(*_fun)(std::vector<float> argv)){
	fun=_fun;
}
/**
 * 初始化种群
 * @param _chrom_num 种群个体数
 * @param _para_num  参数个数
 */
void PSO::crtbp(const int &_chrom_num,const int &_para_num){
	if(_chrom_num>0)chrom_num=_chrom_num;
	if(_para_num>0)para_num=_para_num;
	v=cv::Mat (cv::Size(para_num,chrom_num),CV_32FC1,cv::Scalar(0));
	Population=cv::Mat(cv::Size(para_num,chrom_num),CV_32FC1,cv::Scalar(0));
	ost=std::vector<float>(chrom_num,0);
	post=std::vector<float>(chrom_num,0);
	cv::RNG rng(time(NULL));
	rng.fill(Population, cv::RNG::UNIFORM,bmin ,bmax);
	rng.fill(v, cv::RNG::UNIFORM,0 , 0.5);
	//初始化 Pbest，Gbest
	Population.copyTo(Pbest);
	for(int i=0,min=0;i<chrom_num;i++){
		std::vector<float> argv;
		for(int j=0;j<para_num;j++)
			argv.push_back(Population.at<float>(i,j));
		int result=fun(argv);
		if(i==0||result<min){
			Gbest=argv;
			min=result;
		}
		ost[i]=(result);
	}
}
/**
 * 更新粒子速度、位置
 * 
 */
std::pair<std::vector<float>, float>  PSO::ranking(){
	std::pair<std::vector<float> , float> best(Gbest,fun(Gbest));
	for(int i=0;i<chrom_num;i++){
		std::vector<float> argv;
		for(int j=0;j<para_num;j++){
			//更新速度 updata speed
			v.at<float>(i,j)=(float)(v.at<float>(i,j)+
				c1*(Gbest[j]-Population.at<float>(i,j))*(float)(rand()%100)/100+
				c2*(Pbest.at<float>(i,j)-Population.at<float>(i,j))*(float)(rand()%100)/100    );
			//更新位置 update site
			Population.at<float>(i,j)=(float)Population.at<float>(i,j)+0.5*(float)(v.at<float>(i,j));
			if((float)Population.at<float>(i,j)>bmax)
				Population.at<float>(i,j)=bmax;
			if((float)Population.at<float>(i,j)<bmin)
				Population.at<float>(i,j)=bmin;
			argv.push_back(Population.at<float>(i,j));
		}
		post[i]=fun(argv);
	}
	return best;
}
/**
 * 筛选粒子
 * @param para 0：最小min，1：最大max
 */
void PSO::update(bool para){
	int gbb=fun(Gbest);
	for(int i=0;i<chrom_num;i++){
		if((ost[i]>post[i])^para){
			for(int j=0;j<para_num;j++)
				Pbest.at<float>(i,j)=(float)Population.at<float>(i,j);
			ost[i]=post[i];
		}
		if((post[i]<gbb)^para){
			for(int j=0;j<para_num;j++)
				Gbest[j]=(float)Population.at<float>(i,j);
			gbb=post[i];
		}
	}
}