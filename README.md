# GeneticAlgorithm
这是用C++写的遗传算法，参考《智能算法 30案例分析 第2版》一书，包含TSP、LQR控制器、结合量子算法、多目标优化、粒子群等，由于原作为matlab程序，综合自己思路通过C++写出来，算是练习和开个大坑


- 通过opencv绘制函数曲线图和坐标图
- 一元最优化目标

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
	通过遗传算法求`x^2*sin(3x*pi)`的最大值，增大初始种群数目可加快迭代，增加种群基因编码长度增大迭代稳定性，变异和交叉较小为好

<img src="demo_picture/demo1.png">
	基本第一、二次迭代就接近最优点了，但有点不稳定

## some Problem
	- 在每次选择最优种群个体时，保护当前最优个体加上赌盘选择法可以加快迭代优化，并且增加稳定。
	- C++中构造函数中不能再使用其它重载的构造函数，会失效。