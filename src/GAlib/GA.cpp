#include "GA.h"

GA::GA(const int _chrom_num,const int _gene_num,const float _p_recombin,
	const float _p_mut,const float min,const float max,const int _para_num){
	if(_chrom_num<=0||_gene_num<=0||
		_p_mut<0||_p_mut>1||
		_p_recombin<0||_p_recombin>1||
		max<=min){
		std::cout<<"set GA parameter error,check it!\n";
		exit(1);
	}
	chrom_num=_chrom_num;//种群个体数
	gene_num=_gene_num;//基因数
	p_recombin=_p_recombin;//交叉概率
	p_mut=_p_mut;//变异概率
	search_min=min;//区间下限
	search_max=max;//区间上限
	para_num=_para_num;//编码单元数
}
void GA::solve(float (*_fun)(std::vector<float> argv),const int &_para_num){
	if(_para_num!=0)para_num=_para_num;
	ost=cv::Mat(cv::Size(chrom_num,para_num+1),CV_32FC1,cv::Scalar(0));
	fun=_fun;
}
//输入 种群数，基因数，最小编码，最大编码
//输出生成随机0-1矩阵  Lind<2^32 
cv::Mat GA::crtbp(const int &Nind,const int &Lind,const int&encodemin,const int &encodemax){
	if(Nind>0)chrom_num=Nind;
	if(Lind>0)gene_num=Lind;
	if(gene_num%para_num!=0){std::cout<<"please set (Lind x para_num)!\n";};
	cv::Mat Population(cv::Size(gene_num,chrom_num),CV_8UC1,cv::Scalar(0));
	// cv::randu(Population, 0, 2);//并不随机。。
	cv::RNG rng(time(NULL));
	rng.fill(Population, cv::RNG::UNIFORM,encodemin,encodemax);//UNIFORM or NORMAL
	return Population;
}
//计算适应度
std::pair<std::vector<float>, float> GA::ranking(void){
	//lambada [this]表明是内部类
	std::pair<std::vector<float> , float> best(std::vector<float>(0.0),-1.0/0.0);
	for(int i=0;i<ost.cols;i++){
		std::vector<float> argv;
		for(int j=0;j<para_num;j++)
			argv.push_back((float)ost.at<float>(j,i));
		float m=fun(argv);
		ost.at<float>(ost.rows-1,i)=m;
		if(m>best.second){
			best.first=argv;
			best.second=m;
		}
	}
	return best;
}
//选择优秀个体 bug集中地
GA& GA::select(cv::Mat &Popula,int _method){
	std::vector<std::pair<int, float> > recode_rank_index;
	for(int i=0;i<ost.cols;i++)
		recode_rank_index.push_back(std::pair<int, float>(i,(float)ost.at<float>(ost.rows-1,i)));
	std::sort(recode_rank_index.begin(), recode_rank_index.end(),
					[&](std::pair<int, float> &a,std::pair<int, float> &b){
						if(a.second==b.second)return false;//??不呢为什么
						if((a.second>b.second)^_method)return true;
						else return false;
					});	//排序
	for(int i=0;i<recode_rank_index.size();i++){//修改数据
		recode_rank_index[i].second=(recode_rank_index.size()-i)*(recode_rank_index.size()-i);
	}
	//0^2+1^2+2^2+...+n^2=n(n+1)(2n+1)/6 数学公式很重要的。。。
	unsigned long int sum=(ost.cols*(ost.cols+1)*(2*ost.cols+1))/6;//使用另外一种排序和赌盘选择试试
	// std::for_each(rank.begin(), rank.end(), [&](float & i){sum=sum+std::exp(i);});
	std::vector<uchar> new_Popula;
	//保留最优个体+赌盘选择
	for(int i=0;i<Popula.cols;i++)
			new_Popula.push_back(Popula.at<uchar>(recode_rank_index[0].first,i));
	for(int i=1;i<ost.cols;i++){
		unsigned long int get_one=rand()%sum;//不能使用 float 除数可惜。。
		int _select=0;
		for(unsigned long int add_sum=0;_select<recode_rank_index.size();_select++)
			if(add_sum<get_one)add_sum+=recode_rank_index[_select].second;
			else break;
		for(int j=0;j<Popula.cols;j++)
			new_Popula.push_back(Popula.at<uchar>(recode_rank_index[_select].first,j));
	}
	int rows=Popula.rows,cols=Popula.cols;
	Popula=cv::Mat(new_Popula);
	// cv::resize(Popula, Popula, cv::Size(cols,rows));//使用有问题
	Popula.reshape(0,rows).copyTo(Popula);
	return *this;
}
//交叉  均匀交叉
GA& GA::recombin(cv::Mat &Popula,const float &opt){
	if(opt>0&&opt<1)p_recombin=opt;
	for(int i=0;i<Popula.rows;i++){
		if(rand()%100<p_recombin*100){
			int j=rand()%Popula.rows;
			for(int k=0;k<Popula.cols;k++){
				if(rand()%3<2){
					Popula.at<uchar>(i,k)^=Popula.at<uchar>(j,k);
					Popula.at<uchar>(j,k)^=Popula.at<uchar>(i,k);
					Popula.at<uchar>(i,k)^=Popula.at<uchar>(j,k);
				}
			}	
		}
	}
	return *this;
}
//变异 因为概率建立在统计上，不使用迭代所有成员 均匀变异
GA& GA::mut(cv::Mat &Popula,float opt){
	if(opt>0&&opt<1)p_mut=opt;
	for(int i=p_mut*Popula.cols*Popula.rows;i>0;i--){
		int m=rand()%(Popula.cols*Popula.rows);
		Popula.at<uchar>(m/Popula.cols,m%Popula.cols)^=1;
	}
	return *this;
}
//二进制转十进制 限定区间范围
void GA::bs2rv(cv::Mat &Popula,float min,float max){
	if(min==0)min=search_min;
	if(max==0)max=search_max;
	unsigned long int part=Popula.cols/para_num;
	unsigned long int Max=(unsigned long int)1<<part;Max--;
	//尝试使用at<float>读取 char类型的矩阵？？？
	for(int i=0;i<Popula.rows;i++){
		unsigned long int *sum=new unsigned long int[para_num]();
		for(long int j=0,m=1;j<para_num*part;j++){
			if((bool)Popula.at<uchar>(i,j))sum[j/part]+=m;
			m*=2;
			if(m>Max)m=1;
		}
		for(int j=0;j<para_num;j++)
			ost.at<float>(j,i)=(min+sum[j]*((double)(max-min)/Max)); 
		free(sum);
	}
}
