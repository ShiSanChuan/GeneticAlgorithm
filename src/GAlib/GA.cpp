#include "GA.h"
void GA::initparameter(const int _chrom_num,const int _gene_num,
		const float _p_recombin,const float _p_mut,
		const float min,const float max){
	if(_chrom_num<=0||_gene_num<=0||
		_p_mut<0||_p_mut>1||
		_p_recombin<0||_p_recombin>1||
		max<=min){
		std::cout<<"set GA parameter error,check it!\n";
		exit(1);
	}
	chrom_num=_chrom_num;
	gene_num=_gene_num;
	p_recombin=_p_recombin;
	p_mut=_p_mut;
	search_min=min;
	search_max=max;
}

GA::GA(std::function<float (float&)> _fun,const int _chrom_num,const int _gene_num,const float _p_recombin,
	const float _p_mut,const float min,const float max){
	initparameter(_chrom_num,_gene_num,_p_recombin,_p_mut,min,max);
	fun=_fun;
}
GA::GA(std::function<float (float&,float&)> _fun,const int _chrom_num,const int _gene_num,const float _p_recombin,
	const float _p_mut,const float min,const float max){
	initparameter(_chrom_num,_gene_num,_p_recombin,_p_mut,min,max);
	fun2=_fun;
	
}
//生成随机0-1矩阵  Lind<2^32
cv::Mat GA::crtbp(int Nind,int Lind){
	if(Nind==0)Nind=chrom_num;
	if(Lind==0)Lind=gene_num;
	cv::Mat Population(cv::Size(Lind,Nind),CV_8UC1,cv::Scalar(0));
	// cv::randu(Population, 0, 2);//并不随机。。
	cv::RNG rng(time(NULL));
	rng.fill(Population, cv::RNG::UNIFORM, 0, 2);//UNIFORM or NORMAL
	return Population;
}
//计算适应度
std::pair<float, float> GA::ranking(std::vector<float> &objV){
	//lambada [this]表明是内部类
	std::pair<float, float> best(0,-1.0/0.0);
	std::for_each(objV.begin(), objV.end(),[this,&best](float &i){
							float m=fun(i);
							if(m>best.second){
								best.first=i;
								best.second=m;
							}   
							i=m;});
	return best;
}
//用来求二元最优值
std::vector< std::pair<float, float> > GA::ranking(std::vector<float> &objV,std::vector<float> &objV2){
	std::pair<float, float> best(0,-1.0/0.0);
	std::pair<float, float> best2(0,-1.0/0.0);
	std::vector<float> tem1(objV.size(),-1.0/0.0);
	std::vector<float> tem2(objV2.size(),-1.0/0.0);
	std::vector<std::pair<float, float> > xyz_best;
	for(int i=0;i<objV.size();i++)
		for(int j=0;j<objV2.size();j++){
			float m=fun2(objV[i],objV2[j]);
			if((m>best.second)){
				best.first=objV[i];
				best.second=m;
				best2.first=objV2[j];
				best2.second=m;
			}
			if(tem1[i]<m)tem1[i]=m;
			if(tem2[j]<m)tem2[j]=m;
		}
	xyz_best.push_back(best);
	xyz_best.push_back(best2);
	objV=tem1;
	objV2=tem2;
	return xyz_best ;
}
//选择优秀个体 bug集中地
void GA::select(cv::Mat &Popula,std::vector<float> &rank){
	std::vector<std::pair<int, float> > recode_rank_index;
	for(int i=0;i<rank.size();i++)//创建含下标的rank数据
		recode_rank_index.push_back(std::pair<int, float>(i,rank[i]));
	std::sort(recode_rank_index.begin(), recode_rank_index.end(),
					[&](std::pair<int, float> &a,std::pair<int, float> &b){
						if(a.second>b.second)return true;
						else return false;
					});	//排序
	for(int i=0;i<recode_rank_index.size();i++){//修改数据
		recode_rank_index[i].second=(recode_rank_index.size()-i)*(recode_rank_index.size()-i);
	}
	//0^2+1^2+2^2+...+n^2=n(n+1)(2n+1)/6 数学公式很重要的。。。
	unsigned long int sum=(rank.size()*(rank.size()+1)*(2*rank.size()+1))/6;//使用另外一种排序和赌盘选择试试
	// std::for_each(rank.begin(), rank.end(), [&](float & i){sum=sum+std::exp(i);});
	std::vector<uchar> new_Popula;
	//保留最优个体+赌盘选择
	for(int i=0;i<Popula.cols;i++)
			new_Popula.push_back(Popula.at<uchar>(recode_rank_index[0].first,i));
	for(int i=1;i<rank.size();i++){
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
	cv::resize(Popula, Popula, cv::Size(cols,rows));
}
//交叉  均匀交叉
void GA::recombin(cv::Mat &Popula,float opt){
	if(opt==0)opt=p_recombin;
	for(int i=0;i<Popula.rows;i++){
		if(rand()%100<opt*100){
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
}
//变异 因为概率建立在统计上，不使用迭代所有成员 均匀变异
void GA::mut(cv::Mat &Popula,float opt){
	if(opt==0)opt=p_mut;
	for(int i=opt*Popula.cols*Popula.rows;i>0;i--){
		int m=rand()%(Popula.cols*Popula.rows);
		Popula.at<uchar>(m/Popula.cols,m%Popula.cols)^=1;
	}
}
//二进制转十进制 限定区间范围
std::vector<float> GA::bs2rv(cv::Mat &Popula,float min,float max){
	if(min==0)min=search_min;
	if(max==0)max=search_max;
	std::vector<float> _ost;
	unsigned long int Max=1;
	Max<<=Popula.cols;Max--;
	for(int i=0;i<Popula.rows;i++){
		unsigned long int sum=0;
		for(long int j=0,m=1;j<Popula.cols;j++,m*=2)
			if((bool)Popula.at<uchar>(i,j))sum+=m;
		_ost.push_back(min+sum*((double)(max-min)/Max)); 
	}
	return _ost;
}
