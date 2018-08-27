#include "GA.h"
double GA_TSP::distance(int indexi,int indexj){
	double distance=.0;
	cv::Mat dis1,dis2;
	address.rowRange(indexi, indexi+1).copyTo(dis1);
	address.rowRange(indexj, indexj+1).copyTo(dis2);
	dis1=dis1-dis2;
	dis1=dis1.mul(dis1);
	for(int i=0;i<dis1.cols;i++)
		distance+=dis1.at<float>(0,i);
	return cv::sqrt(distance);
}
//传入目的解
void GA_TSP::TSPsolve(cv::Mat &_address){
	address=_address;
	gene_num=para_num=_address.rows;//因为基因的编码没有改变
	ost=cv::Mat(cv::Size(chrom_num,1),CV_32FC1,cv::Scalar(0));
}
//对成员进行评分
std::pair<std::vector<float>, float> GA_TSP::ranking(cv::Mat &_Popula){
	std::pair<std::vector<float> , float> best(std::vector<float>(0.0),-1.0/0.0);
	for(int i=0;i<_Popula.rows;i++){
		double _distance=0.0;
		for(int j=1;j<_Popula.cols;j++)
			_distance+=distance(_Popula.at<uchar>(i,j-1),_Popula.at<uchar>(i,j));
		_distance+=distance(_Popula.at<uchar>(i,0),_Popula.at<uchar>(i,_Popula.cols-1));
		ost.at<float>(0,i)=1/_distance;
		if(1/_distance>best.second){
			best.second=1/_distance;
			std::vector<float> argv;
			for(int j=0;j<_Popula.cols;j++)
				argv.push_back((uchar)_Popula.at<uchar>(i,j));
			best.first=argv;
		}
	}
	return best;
}
//创建种群
cv::Mat GA_TSP::crtbp(int encodemax){
	if(encodemax>0)gene_num=para_num=encodemax;
	cv::Mat Population(cv::Size(1,chrom_num),CV_32FC1,cv::Scalar(1));
	cv::Mat gene(cv::Size(gene_num,1),CV_32FC1,cv::Scalar(0));
	for(int i=0;i<gene.cols;i++)
			gene.at<float>(0,i)=i;
	Population*=gene;//创建一个标准种群
	Population.convertTo(Population, CV_8UC1);
	return Population;
}
//反向旋转 
GA_TSP& GA_TSP::recombin(cv::Mat &Popula,const float &opt){
	if(opt>0&&opt<1)p_recombin=opt;
	uchar * tem=new uchar[2*Popula.cols];
	for(int i=1;i<10;i++){
		uchar *gene1=Popula.ptr<uchar>(0);
		uchar *gene2=Popula.ptr<uchar>(Popula.rows-i);
		memcpy(gene2,gene1,sizeof(uchar)*Popula.cols);
	}

	for(int i=0;i<Popula.rows;i++){
		if(rand()%100<p_recombin*100){
			uchar *gene=Popula.ptr<uchar>(i);
	 		memcpy(tem,gene,sizeof(uchar)*Popula.cols);
	 		memcpy(tem+Popula.cols,gene,sizeof(uchar)*Popula.cols);
			memcpy(gene,tem+rand()%Popula.cols,sizeof(uchar)*Popula.cols);
		}
	}
	delete  []tem;
	return *this;
}
//变异
GA_TSP& GA_TSP::mut(cv::Mat &Popula,float opt){
	if(opt>0&&opt<1)p_mut=opt;
	for(int i=1;i<Popula.rows;i++){
		if(rand()%100>opt*100)
		for(int j=0;j<p_mut*Popula.cols;j++){
			uchar *p=Popula.ptr<uchar>(i);
			uchar *a=p+rand()%Popula.cols;
			uchar *b=p+rand()%Popula.cols;
			while((*b)==(*a))b=p+rand()%Popula.cols;
			(*a)=(*a)^(*b);
			(*b)=(*a)^(*b);
			(*a)=(*a)^(*b);
		}
	}
	return *this;
}