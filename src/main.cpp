#include <iostream>
#include <algorithm>
#include <vector>
// #include <iomanip> //控制cout 输出格式
#include "cvplot.h"
#include "opencv2/opencv.hpp"
#include "GA.h"
//Y=sin(X)
int main(int argc, const char** argv)
{
	srand(time(NULL));
	cv::Mat Popula;
	std::vector<float> ost;
	std::vector<float> data;
	std::vector<std::pair<float, float> > recode;
	GA ga([](float &i)->float{return (i*i*std::sin(i*3*pi));});
	//种群越大，越稳定，基因编码越多迭代越稳定
	Popula=ga.crtbp();
	for(int i=0;i<20;i++){
		ost=ga.bs2rv(Popula);

		recode.push_back(ga.ranking(ost));
		data.push_back(recode[recode.size()-1].second);
		ga.select(Popula, ost);
		ga.recombin(Popula);
		ga.mut(Popula);
	}
	//绘图
	{
		auto name="math";
		cvplot::setWindowTitle(name,"origin curve");
		cvplot::moveWindow(name, 0, 0);
		cvplot::resizeWindow(name, 400, 400);
		auto &figure=cvplot::figure(name);
		std::vector<std::pair<float, float> > data;
		for(float i=0;i<10;i+=0.025)
			data.push_back({i,i*i*std::sin(i*3*pi)});
		figure.series("x^2*sin(3x*pi)").set(data).color(cvplot::Green);
		figure.series("recode").set(recode).type(cvplot::Dots).color(cvplot::Red);
		figure.border(30).show(false);
	}
	{
		auto name="GA";
		cvplot::setWindowTitle(name,"GA");
		cvplot::moveWindow(name, 400, 0);
		cvplot::resizeWindow(name, 400, 400);
		auto &figure=cvplot::figure(name);
		figure.series("count").setValue(data).color(cvplot::Orange);
		figure.border(30).show();
	}
	cv::waitKey(0);
	return 0;
}
