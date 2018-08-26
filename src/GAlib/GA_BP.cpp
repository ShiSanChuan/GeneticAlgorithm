#include "GA.h"
//y=fun(x*w+b)
void GA_BP::BPsolve(cv::Mat &_input,cv::Mat &_output){
	input=_input;
	output=_output;
	implication_num=input.cols*2+1;
	para_num=implication_num*input.cols+implication_num;
	para_num+=implication_num*output.cols+output.cols;
	ost=cv::Mat(cv::Size(chrom_num,para_num+1),
				CV_32FC1,cv::Scalar(0));
}
std::pair<std::vector<float>, float> GA_BP::ranking(void){
	std::pair<std::vector<float> , float> best(std::vector<float>(0.0),1.0/0.0);
	cv::Mat in_imp,in_imp_b,imp_out,imp_out_b,err;
	for(int i=0;i<ost.cols;i++){
		ost.colRange(i,i+1).rowRange(0,implication_num*input.cols).copyTo(in_imp);
		ost.colRange(i,i+1).rowRange(implication_num*input.cols,
										implication_num*(input.cols+1)).copyTo(in_imp_b);
		ost.colRange(i,i+1).rowRange(implication_num*(input.cols+1),
				implication_num*(input.cols+1)+implication_num*output.cols).copyTo(imp_out);
		ost.colRange(i,i+1).rowRange(
			implication_num*(input.cols+1)+implication_num*output.cols,
			para_num).copyTo(imp_out_b);
		float _err=0;
		for(int j=0;j<input.rows;j++){
			err=(input.rowRange(j, j+1)*in_imp.reshape(0,input.cols)+in_imp_b.reshape(0,1))*
									imp_out.reshape(0,implication_num)+imp_out_b.reshape(0,1)-
									output.rowRange(j,j+1);
			cv::multiply(err, err, err);
			for(int k=0;k<err.cols;k++)
				_err+=(float)err.at<float>(0,k);
		}
		ost.at<float>(ost.rows-1,i)=std::sqrt(_err);
		if(_err<best.second)best.second=_err;
	}
	std::cout<<ost.rowRange(ost.rows-1, ost.rows)<<std::endl;
	return best;
}