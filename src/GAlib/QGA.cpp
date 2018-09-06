#include "GA.h"
//量子编码的遗传算法
//由于原编码为 CV_8UC1 转为 CV_32FC1 行内缩短4倍

/**
 * 创建量子种群
 * @param Nind      列数
 * @param Lind      行数
 * @param encodemin 最小编码
 * @param encodemax 最大编码
 */
cv::Mat QGA::crtbp(const int &Nind,const int &Lind){
	if(Nind>0)chrom_num=Nind;
	if(Lind>0)gene_num=Lind;
	if(gene_num%para_num!=0){std::cout<<"please set (Lind x para_num)!\n";};
	cv::Mat Population(cv::Size(gene_num,chrom_num),CV_8UC1,cv::Scalar(0));
	// cv::randu(Population, 0, 2);//并不随机QGA。。
	for(int i=0;i<Population.rows;i++)
		for(int j=0;j<Population.cols/4;j++)
			Population.at<float>(i,j)=pi/4;//初始化都为45'
	return Population;
}

/**
 * 将数据编码
 * @param Popula [种群基因]
 * @param min    [区间下限]
 * @param max    [区间上限]
 */
void QGA::bs2rv(cv::Mat &Popula,float min,float max){
	if(min==0)min=search_min;
	if(max==0)max=search_max;
	unsigned long int part=Popula.cols/para_num;//每个变量基因位数
	part/=4;//char 8 -> float 32
	unsigned long int Max=(unsigned long int)1<<part;Max--;
	//尝试使用at<float>读取 char类型的矩阵？？？可行
	for(int i=0;i<Popula.rows;i++){
		unsigned long int *sum=new unsigned long int[para_num]();
		for(long int j=0,m=1;j<para_num*part;j++){
			double a=std::sin((float)(Popula.at<float>(i,j)));
			if(rand()%100>a*a*100)sum[j/part]+=m;//原书P82描述为平方项
			m*=2;
			if(m>Max)m=1;
		}
		for(int j=0;j<para_num;j++)
			ost.at<float>(j,i)=(min+sum[j]*((double)(max-min)/Max)); 
		free(sum);
	}
	// std::cout<<ost<<"\t"<<ost.cols<<"\t"<<ost.rows<<std::endl;
}
/**
 * select 通过量子门旋转更新群体 更新方法: P82
 *x	best 	f(x)>f(best)	alpha	ab>0	ab<0	a=0	b=0
 *0	0		false			0		0		0		0	0
 *0	0		true			0		0		0		0	0
 *0	1		false			0.01pi	1		-1		0	+-1
 *0	1		true			0.01pi	-1		1		+-1	0
 *1	0		false			0.01pi	-1		1		+-1	0
 *1	0		true			0.01pi	1		-1		0	+-1
 *1	1		false			0		0		0		0	0
 *1	1		true			0		0		0		0	0
 *其实不用这么复杂，即角度想最优个体调整
 *ab>0时，0~pi/2 || -pi~-pi/2  ab<0: -pi/2~0 || pi/2~pi
 */
QGA& QGA::select(cv::Mat &Popula){
	int best_index=0;
	int sign=0;//旋转方向
	double delta=0.01*pi;//旋转角度
	float best=(float)ost.at<float>(ost.rows-1,0);
	for(int i=1;i<ost.cols;i++)
		if((float)ost.at<float>(ost.rows-1,i)>best){
			best_index=i;
			best=(float)ost.at<float>(ost.rows-1,i);
		}
	for(int i=0;i<Popula.rows;i++)
		for(int j=0;j<Popula.cols/4;j++){
			double a=std::sin((float)(Popula.at<float>(i,j)));
			double b=std::sin((float)(Popula.at<float>(best_index,j)));
			if((rand()%100>a*a*100)^(rand()%100>b*b*100)){
				float angle=(float)(Popula.at<float>(i,j));
				if((float)(ost.at<float>(ost.rows-1,i))>(float)(ost.at<float>(ost.rows-1,best_index))){
					// if()
					if((angle>0&&angle<pi/2)||(angle>-pi&&angle<-pi/2))
						sign=-1;
					else sign=1;
				}else{
					if((angle>0&&angle<pi/2)||(angle>-pi&&angle<-pi/2))
						sign=1;
					else sign=-1;
				}
				Popula.at<float>(i,j)=angle+sign*delta;
			}else delta=0;
		}

	return *this;
}