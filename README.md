# GeneticAlgorithm
这是用C++写的遗传算法，参考《智能算法 30案例分析 第2版》一书，包含TSP、LQR控制器、结合量子算法、多目标优化、粒子群等，由于原作为matlab程序，综合自己思路通过C++写出来，算是练习和开个大坑


- 通过opencv绘制函数曲线图和坐标图
- 一元最优化目标
- 二元函数优化目标
- 基于遗传算法的BP神经网络（施工中）

## How to use
```
git clone https://github.com/ShiSanChuan/GeneticAlgorithm.git
cd GeneticAlgorithm/
cmake .
make -j4
cd src/
./GA
```

## Recode
### 一元函数优化
- 通过遗传算法求`x^2*sin(3x*pi)`的最大值，增大初始种群数目可加快迭代，增加种群基因编码长度增大迭代稳定性，变异和交叉较小为好;
- 对于遗传算法中的赌盘轮巡法，最常见直接计算所有个体函数的累加值作为随机值的最大值，但因为数据中可能有负数，所以将所有数据减去这个最小值，但这样结果会使中间的数据频繁出现，无法很好的表现最优值;

<img src="demo_picture/demo1_1.png">

- 若再将fun(x)后处理的数据再次进行exp(x),可以去除有负数累加的问题，但若目标函数的最优值与次最优值很近，或者fun(x)的数据本身集中在[-inf,1]，也无法很好的区分最优值;

<img src="demo_picture/demo1_2.png">

- 若将fun(x)后的数据进行排序，在将排序的大小换为整数1,...,n，在通过指数exp或者平方处理，可以很好区分最优值。

<img src="demo_picture/demo1_3.png">

- 因为种群编码长度受计算机影响（64位），因此搜索区间太大会使精度下降，因此进行一次遗传算法寻求全局最优后再次缩小范围求获得的精度更高。

### 二元函数优化
	与一元函数优化基本类似，不过在rank中需要在二元中需找对应最大解。

<img src="demo_picture/demo2_1.png">
<img src="demo_picture/demo2_2.png">


## some Problem
	- 在每次选择最优种群个体时，保护当前最优个体加上赌盘选择法可以加快迭代优化，并且增加稳定。
	- C++中构造函数中不能再使用其它重载的构造函数，会失效。
	- cv::Mat 的结构体，使用at时，返回的是一个template的量，因此在使用逻辑表达式的时候最好先强制转换一下。
	- 使用cstdarg，C++可变参数函数：

```cpp
#include <cstdarg>
void testarg(count,...){
	va_list ap;
	va_start(ap, count);
	for(int i=0;i<count;i++)
	std::cout<<(int)va_arg(ap, int)<<std::endl;;
	va_end(ap);
	return;
}
```