#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#include<io.h>
#include<direct.h>    //头文件  

#include "vrp_5_2.h"
#include<vector>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<string>
using namespace std;

/*
参数列表说明：
vector<double> a		:	是一个规一化的无序数组，和为1
返回值：
返回数组下标索引
*/
int RouletteWheelSelection(vector<double> a)
{
	double fSlice = rand() % 1000 / (double)1001;

	double cfTotal = 0.0; //概率统计

	int SelectedGenomePos = 0; //被选中的个体下标
	int i;
	for (i = 0; i < a.size(); ++i)
	{

		cfTotal += a[i];
		if (cfTotal > fSlice)
		{
			SelectedGenomePos = i;
			break;
		}
	}

	return SelectedGenomePos;
}

/*
功能:初始化客户点、客户点之间的行驶时间和行驶距离
distance[]:距离的文件路径字符数组
time[]:时间的文件路径字符数组
spec[]:
算法步骤：
1.通过文件的路径，建立距离文件、时间文件和种群文件的三个文件流对象
2.从r_spec文件流中读取数据，为250个客户点赋值
3.关闭r_spec文件流
4.从r_time、r_distance文件流中读取数据，为客户点i与客户点j的行驶时间、行驶距离赋值
5.为客户点251（表车场）到其它客户点[0,250]的行驶时间、行驶距离赋值，因为客户点0与客户点n+1都是车场
6.关闭时间流和距离流
*/
void getData(char distance[], char time[], char spec[]) {

	//打开文件读取文件内容到文件流对象 in中
	//1.建三个文件流对象，和文件关联
	fstream r_distance(distance, ios::in), r_time(time, ios::in), r_spec(spec, ios::in);
	//cust_num=50;

	char line[300];
	for (int i = 0; i < 4; ++i)
	{
		r_spec.getline(line, 300);			//把文件流中的300个字节的数据放入line中。
		cout << line << endl;
	}

	for (int i = 0; i<2; ++i)
	{
		r_spec >> max_capacity;				//把文件流中的两个int类型给max_capacity
		cout << max_capacity << endl;
	}

	for (int i = 0; i < 4; ++i)
	{
		r_spec.getline(line, 300);
		cout << line << endl;
	}
	//2.从r_spec文件流中读取数据为250个客户点赋值
	for (int i = 0; i <= cust_num; ++i) {			//为客户点0-249，赋值
		r_spec >> customer[i].id;				//客户点i的序号
		r_spec >> customer[i].x;				//客户点i的x值
		r_spec >> customer[i].y;				//客户点i的y值
		r_spec >> customer[i].demand;			//客户点i的需求量
		r_spec >> customer[i].b_time;			//客户点i的最早服务时间
		r_spec >> customer[i].e_time;			//客户点i的最迟服务时间
		r_spec >> customer[i].s_time;			//客户点i的服务时间
		cout << " id:" << customer[i].id << " x:" << customer[i].x << " y:" << customer[i].y << " demand:" << customer[i].demand << " b_time:" << customer[i].b_time << "  e_time:" << customer[i].e_time << "  s_time:" << customer[i].s_time << endl;
	}
	//3.关闭r_spec文件流
	r_spec.close();

	//4.从r_time、r_distance文件流中读取数据，为客户点i与客户点j的行驶时间、行驶距离赋值
	for (int i = 0; i <= cust_num; ++i)
		for (int j = 0; j <= cust_num; ++j) {
			r_time >> peer_time[i][j];			//车辆从客户点i到客户点j的行驶时间
			r_distance >> peer_distance[i][j];	//车辆从客户点i 到客户点j的行驶距离
		}

	//5.为客户点0-250与客户点251的行驶时间与行驶距离赋值，因为客户点0与客户点251都表示车场
	for (int i = 0; i <= cust_num; ++i) {				//为客户点251（即：车场 0）添加客户点i到客户点cust_num + 1 = 251的行驶时间、行驶距离
		peer_time[i][cust_num + 1] = peer_time[i][0];
		peer_time[cust_num + 1][i] = peer_time[0][i];
		peer_distance[i][cust_num + 1] = peer_distance[i][0];
		peer_distance[cust_num + 1][i] = peer_distance[0][i];
	}
	customer[cust_num + 1] = customer[0];			//索引为cust_num + 1与 0的都表示车场
													//6.关闭时间流和距离流
	r_time.close();								//时间流关闭
	r_distance.close();							//距离流关闭
}


//清空路径：容量、延迟时间、行驶距离、行驶时间、等待时间；清空路径中所有的客户点
void clear_route(Route &route) {
	route.capacity = 0;
	route.delay_time = 0;
	route.travel_dist = 0;
	route.travel_time = 0;
	route.wait_time = 0;
	route.customers.clear();
}

/*清空解：
1.清空解中的路径；
2.将一个解中的5个目标的值置为：0；
*/
void clean_chromosome(Chromosome &chrome) {

	chrome.routes.clear();						//清空解中的路径
	for (int i = 0; i<FUNC_NUM; ++i)
		chrome.f[i] = 0;							//将一个解中的5个目标的值置为：0

}

//计算解的目标函数：f[0]:车辆数目；f[1]:总行驶距离；f[2]:最长路径行驶时间；f[3]:总等待时间；f[4]:总延迟时间
void compute_f(Chromosome &chrome) {

	for (int i = 0; i<FUNC_NUM; ++i)
		chrome.f[i] = 0;
	chrome.f[0] = chrome.routes.size();
	for (int i = 0; i<chrome.routes.size(); ++i)
		chrome.f[1] += chrome.routes[i].travel_dist;
	int index = 0; double travel_time = 0;
	for (int i = 0; i<chrome.routes.size(); ++i)
		if (chrome.routes[i].travel_time>travel_time)
		{
			travel_time = chrome.routes[i].travel_time;
			index = i;
		}
	chrome.f[2] = chrome.routes[index].travel_time;
	chrome.cub_len = 0;

	//  chrome.f[1]=0;
	//chrome.f[2]=0;
	for (int i = 0; i<chrome.routes.size(); ++i) chrome.f[3] += chrome.routes[i].wait_time;
	for (int i = 0; i<chrome.routes.size(); ++i) chrome.f[4] += chrome.routes[i].delay_time;

}

//计算路径的长度，b_distance[*iter]、a_distance[*iter]的计算，为LS6提供数据
void compute_route_distance(Route &route) {                                               //compute the traval distance for a route

	route.travel_dist = 0;
	vector<int>::iterator iter;
	for (iter = route.customers.begin(); iter != route.customers.end() - 1; ++iter) {
		route.travel_dist += peer_distance[*iter][*(iter + 1)];
		a_distance[*iter] = 0;
		b_distance[*iter] = 0;
	}
	iter = route.customers.begin();
	b_distance[*iter] = 0;
	a_distance[*iter] = 0;

	//i之前的计算
	for (iter = route.customers.begin() + 1; iter != route.customers.end() - 1; ++iter)
		b_distance[*iter] += (b_distance[*(iter - 1)] + peer_distance[*(iter - 1)][*iter]);
	//i之后的计算
	for (iter = route.customers.end() - 2; iter >= route.customers.begin() + 1; --iter)
		a_distance[*iter] += (a_distance[*(iter + 1)] + peer_distance[*iter][*(iter + 1)]);

}

//计算路径中各个客户点的行驶时间、等待时间、延迟时间和最大合法延迟时间（方便插入客户点）
void compute_route_time(Route &route) {                                                   //compute the traval time, wait time and delay time for a route
	route.travel_time = 0;
	route.delay_time = 0;
	route.wait_time = 0;
	vector<int>::iterator iter;
	total_wait[0] = 0;
	total_delay[0] = 0;
	for (iter = route.customers.begin(); iter != route.customers.end() - 1; ++iter) {
		double wait = 0, delay = 0;
		route.travel_time += peer_time[*iter][*(iter + 1)];
		a_time[*(iter + 1)] = route.travel_time;
		//求延迟时间
		if (route.travel_time>customer[*(iter + 1)].e_time)	delay = route.travel_time - customer[*(iter + 1)].e_time;			//求客户点（1--- n - 1）延迟时间
		route.delay_time += delay;
		total_delay[*(iter + 1)] = total_delay[*iter] + delay;								//总的下一个点的延迟时间 = 本客户点的延迟时间 + 下一个点的延迟时间
																							//求等待时间
		if (route.travel_time<customer[*(iter + 1)].b_time)	wait = customer[*(iter + 1)].b_time - route.travel_time;
		w_time[*(iter + 1)] = wait;
		total_wait[*(iter + 1)] = total_wait[*iter] + wait;
		route.wait_time += wait;
		route.travel_time += (wait + customer[*(iter + 1)].s_time);							//路径的总行驶时间 = [k + 1]的行驶时间 + [k + 1]的等待时间 + [k + 1]的服务时间；
		l_time[*(iter + 1)] = route.travel_time;											//第 k + 1的离开时间 = 第k + 1的行驶时间
	}
	iter = route.customers.begin();
	l_time[*iter] = 0;
	//求各个客户点的最大合法延迟时间
	max_wait[*iter] = customer[*iter].e_time - a_time[*iter];								//车场最大合法延迟时间 = 车场关闭时间 - 到达车场的时间

	for (iter = route.customers.end() - 2; iter != route.customers.begin(); --iter) {			//对于k......1的客户点的最大合法延迟时间
																								//迟到																		早到
		max_wait[*iter] = min(customer[*iter].e_time + max_delay_time - a_time[*iter], w_time[*iter] + max_wait[*(iter + 1)]);
	}
	iter = route.customers.begin();
	total_wait[*iter] = 0;
	total_delay[*iter] = 0;

}

//计算路径route的容量
void compute_route_capacity(Route &route) {                                                      //compute the route capacity

	route.capacity = 0;
	vector<int>::iterator iter;
	for (iter = route.customers.begin() + 1; iter != route.customers.end() - 1; ++iter) {
		route.capacity += customer[*iter].demand;
		b_demand[*iter] = 0;
		a_demand[*iter] = 0;
	}
	iter = route.customers.begin();
	a_demand[*iter] = 0;
	b_demand[*iter] = 0;
	//计算路径中各个客户点i在路径之前的容量b_demand,和之后的容量a_demand(包含i的需求)
	for (iter = route.customers.begin() + 1; iter != route.customers.end() - 1; ++iter)
		b_demand[*iter] += (b_demand[*(iter - 1)] + customer[*iter].demand);
	for (iter = route.customers.end() - 2; iter >= route.customers.begin() + 1; --iter)
		a_demand[*iter] += (a_demand[*(iter + 1)] + customer[*iter].demand);

}

/*
初始化路径,产生新解：
分别将请求数组中的元素插入到路径中，插入完，产生新解；此处没有求:delay_time
由于刚开始的客户点的选择（通过234行：int r_ind=rand()%cust.size();）产生，所有生成的解是随机的
*/
void init_route(Chromosome &chrome) {
	clean_chromosome(chrome);
	vector<int> cust, temp_cust;						//声明:请求所在的集合
														//生成服务请求集合
	for (int i = 1; i <= cust_num; ++i)					//[0,250]
		cust.push_back(i);							//该函数将一个新的元素加到vector的最后面

	Route new_route;								//生成一个路径
	while (true) {
		clear_route(new_route);						//清空新的路径 
		new_route.customers.push_back(0);			//把车场编号0添加到当前路径当中
		double w_time = 0;							//临时变量等待时间初始化为0					
		double t_time = 0;							//临时变量行驶时间初始化为0	
		double capacity = 0;							//临时变量容量初始化为0	
		int last = 0;									//记录上一个客户点编号							
		while (true) {                                                                                                         //create a route
			vector<int>::iterator iter = cust.begin();
			int r_ind = rand() % cust.size();			//随机选择一个客户点的请求,由于是随机选择客户点，所有产生的解是随机的
			t_time += peer_time[last][cust[r_ind]];	//路径时间等于路径上行驶的时间
			w_time = t_time > customer[cust[r_ind]].b_time ? 0 : customer[cust[r_ind]].b_time - t_time;								//根据公式求等待时间
																																	//根据容量约束或时间窗约束来判断，是否可以加入当前路径的情况:
																																	//不可以加入当前路径
			if (capacity + customer[cust[r_ind]].demand > max_capacity ||		//加入之后，容量大于路径最大容量；
				t_time > customer[cust[r_ind]].e_time + max_delay_time ||	//迟到，行驶时间大于客户点r_ind的最迟时间 + 最大延迟时间
																			//将客户点插入到第i个点时到车场0的时间,即最后到达车场的时间大于车场的最后工作时间
				t_time + w_time + customer[cust[r_ind]].s_time + peer_time[cust[r_ind]][0] > customer[0].e_time) {

				t_time -= peer_time[last][cust[r_ind]];//将行驶时间归于上一个行驶时间
				temp_cust.push_back(cust[r_ind]);	//将这个客户点请求放入临时客户点集
				cust.erase(iter + r_ind);				//不能用的客户点从cust客户点集中删除
				if (cust.size() == 0)
					break;			//客户点数为0时,跳出循环
				continue;							//否则,继续判断
			}
			else {		//满足时间窗约束与车辆容量约束的情况:
				t_time += w_time + customer[cust[r_ind]].s_time;		//行驶时间 = 上一个行驶时间 + 在r_ind 客户点的等待时间 + 在r_ind客户点的服务时间
				capacity += customer[cust[r_ind]].demand;				//第r_ind客户点的容量 = 前一个客户点的容量 + 第r_ind客户点的容量
				new_route.customers.push_back(cust[r_ind]);			//将客户点插入到路径中
																	//此处可以求延迟时间
				last = cust[r_ind];
				cust.erase(iter + r_ind);								//从集合cust中删除第r_ind个客户点
				for (iter = cust.begin(); iter != cust.end(); ++iter)		//添加成功后，将剩余的客户点加入到temp_cust客户点集中
					temp_cust.push_back(*iter);						//将cust集合中的元素copy到temp_cust集合中
				cust = temp_cust;										//将剩余的服务请求重新加入到cust集合中					
				temp_cust.clear();									//临时集合清空

			}

			if (cust.size() == 0)									//客户点数为0时,跳出循环;当cust集合中的元素为0是因为cust中的客户点都不适合添加到此路径,可以开始新路径的创建
				break;
		}

		cust = temp_cust;
		temp_cust.clear();
		new_route.customers.push_back(0);							//路径的最后一个点0的插入
		compute_route_distance(new_route);							//计算新路径的行驶距离
		compute_route_time(new_route);
		compute_route_capacity(new_route);
		chrome.routes.push_back(new_route);
		if (cust.size() == 0)
			break;
	}
}

/*将存档中的解与新解进行比较
chrome_1:
chrome_2:
chrome_2的关联矩阵的值都大于chrome_1，better = true
只要chrome_2的关联矩阵有一个方向比chrome_1小，2好，1不好，返回false
*/
bool is_better(Chromosome &chrome_1, Chromosome &chrome_2) {

	bool better = false;
	for (int i = 0; i<FUNC_NUM; ++i) {
		if (chrome_2.f[i] - chrome_1.f[i]>0.000001)
			better = true;
		if (chrome_1.f[i] - chrome_2.f[i]>0.000001)
			return false;	//只要chrome_1的关联矩阵有一个方向比旧解大，chrome_1是坏的，chrome_2有可取之处，返回false
	}//chrome_2的关联矩阵的值都大于chrome_1，better = true，chrome_1是好的，chrome_2是坏的
	return better;
}

/*
功能：比较chrome_1与chrome_2的优劣
所有方向上chrome_1比chrome_2时，chrome_1才是好的
只要chrome_1有一个方向上比chrome_2差，就返回false,即：chrome_1不优于chrome_2
*/
bool is_better_box(Chromosome &chrome_1, Chromosome &chrome_2) {

	bool better = false;
	for (int i = 0; i<FUNC_NUM; ++i) {
		if (chrome_2.box_f[i] - chrome_1.box_f[i]>0.000001)				//所有方向上chrome_1比chrome_2时，chrome_1才是好的				
			better = true;
		if (chrome_1.box_f[i] - chrome_2.box_f[i]>0.000001)				//只要chrome_1有一个方向上比chrome_2差，就返回false,即：chrome_1不优于chrome_2
			return false;
	}
	return better;

}


/*
功能：用LS6,检查最优的插入位置,没有进行两路径序列的交换：
返回true:找到最优的位置了
返回false:没找到最优位置
check_best_pos(chrome.routes[r_ind],c_ind,chrome.routes[i],c_pos,delta,r_lamda,chrome.routes.size(),dist,now_time,now_wait,now_delay)
check_best_pos(Route &org_route, int cust_pos, Route &route, int &c_pos, double &min_delta, double r_lamda[], int r_num, double now_dist,double now_time,double now_wait, double now_delay)
org_route:	chrome.routes[r_ind]，解chrome的第r_ind条路径
cust_pos：	c_ind，随机从[0,r_ind.size - 2]中产生的索引
route：		chrome.routes[i]，解chrome的第i条路径,遍历的外层循环
c_pos：		c_pos			在函数中，只赋值，用于记录最优插入位置
min_delta：	delta = 0		变量增量
r_lamda[]：	r_lamda
r_num：		chrome.routes.size()，此解中路径的数目
now_dist：	dist，除org_route(r_ind)和route(i) 路径的路径长度之和
now_time:	now_time 	//记录一条路径的行驶时间,当优化f[2]、f[4]、f[5]，now_time = 0;当优化f[3]，now_time = 总行驶时间最长的路径的行驶时间
now_wait：	now_wait,除org_route(r_ind)和route(i) 路径的等待时间之和
now_delay：	now_delay,除org_route(r_ind)和route(i) 路径的延迟时间之和

结果：
用min_dalta记录LS6领域操作后对应目标的值，其中对于f[2]，min_delta记录的是r_ind和i路径中最长路径行驶时间
即：找到最优插入位置后，对于f[0]是最少车辆数目；
对于f[1]是总行驶距离；
对于f[2]是r_ind路径与i路径的最长行驶时间
对于f[3]是总等待时间
对于f[4]是总延迟时间
用c_pos记录最优的插入位置的索引
feasible记录是否找到插入位置
*/
bool check_best_pos(Route &org_route, int cust_pos, Route &route, int &c_pos, double &min_delta, double r_lamda[], int r_num, double now_dist, double now_time, double now_wait, double now_delay) {    //2-opt*

	bool feasible = false;
	min_delta = 9999999;												//设置成一个大值
																		//LS6操作后的：
																		//i路径的距离增量,
	double	delta_dist = 0;
	//r_ind路径的时间增量,
	double	delta_time_org = 0;
	//i路径的时间增量，
	double	delta_time = 0;
	//增量
	double	delta = 0;
	//r_ind路径距离，	
	double	dist_org;
	//i路径距离，	
	double	dist_mod;
	//r_ind路径时间，
	double	time_org;
	//i路径时间
	double	time_mod;

	Route temp_route = route;																					//chrome.routes[i]，解chrome的第i条路径,遍历的外层循环
	temp_route.customers[0] = cust_num + 1;																	//251,客户点的标号
	temp_route.customers[temp_route.customers.size() - 1] = cust_num + 1;										//251
	compute_route_capacity(temp_route);
	compute_route_time(temp_route);
	compute_route_distance(temp_route);
	/*
	容量约束检查：
	r_ind路径中的cust_pos客户点						 i路径的客户点客户点迭代器
	*/
	vector<int>::iterator	org_c_iter = org_route.customers.begin() + cust_pos, c_iter;
	//将r_ind路径cust_pos客户点后的序列与i路径适当c_iter客户点后的序列进行交换
	for (c_iter = temp_route.customers.begin() + 1; c_iter != temp_route.customers.end(); ++c_iter) {
		//两种情况不行，指针向下一位置移动
		//对于r_ind路径：if cust_ind客户点之前容量 + c_iter之后容量 > max_capacity，++c_iter
		if (b_demand[*org_c_iter] + a_demand[*c_iter] > max_capacity)
			continue;
		//否则，r_ind路径交换容量检查通过，检查i路径
		//对于i路径，路径i只剩c_iter - 1之前部分，r_ind路径只剩cust_ind + 1之后部分，连接检查容量约束
		if (b_demand[*(c_iter - 1)] + a_demand[*(org_c_iter + 1)] > max_capacity)
			continue;
		/*对于i路径连接后，容量检查通过，进行时间检查
		第1种时间检查：看客户点是不是0
		超时两种情况，即：c_iter或cust_ind +　１	没有到0
		到0
		*/
		if (*c_iter != 0) {
			if (l_time[*org_c_iter] + peer_time[*org_c_iter][*c_iter]>customer[*c_iter].e_time + max_delay_time)
				continue;
		}	//*c_iter == 0
		else {
			if (l_time[*org_c_iter] + peer_time[*org_c_iter][*c_iter]>customer[*c_iter].e_time) continue;
		}

		if (*(org_c_iter + 1) != 0) {  //路径r_ind客户点c_ind的下一点不是0
			if (l_time[*(c_iter - 1)] + peer_time[*(c_iter - 1)][*(org_c_iter + 1)]>customer[*(org_c_iter + 1)].e_time + max_delay_time)continue;
		}
		else {
			if (l_time[*(c_iter - 1)] + peer_time[*(c_iter - 1)][*(org_c_iter + 1)]>customer[*(org_c_iter + 1)].e_time)continue;
		}
		/*
		第二种时间检查，看i路径的c_iter后序列插入r_ind路径cust_ind后   和   将r_ind路径cust_ind + 1后序列插入i路径c_iter - 1前
		*/
		//将c_iter后序列插入r_ind路径，到c_iter时的，时间增量
		delta_time_org = l_time[*org_c_iter] + peer_time[*org_c_iter][*c_iter] - a_time[*c_iter];
		if (delta_time_org > max_wait[*c_iter])	//不能插入
			continue;
		//c_iter后序列插入r_ind路径可以插入
		//将cust_ind + 1后序列插入i路径，到org_c_iter+1时的，时间增量
		delta_time = l_time[*(c_iter - 1)] + peer_time[*(c_iter - 1)][*(org_c_iter + 1)] - a_time[*(org_c_iter + 1)];
		if (delta_time>max_wait[*(org_c_iter + 1)])
			continue;
		//两种时间检查后
		//求交换后，r_ind路径距离之和dist_org
		dist_org = b_distance[*org_c_iter] + peer_distance[*org_c_iter][*c_iter] + a_distance[*c_iter];
		//求交换后，i路径距离之和dist_mod
		dist_mod = b_distance[*(c_iter - 1)] + peer_distance[*(c_iter - 1)][*(org_c_iter + 1)] + a_distance[*(org_c_iter + 1)];

		if (dist_org == 0) {
			--r_num;					//解的路径数减1
		}
		if (dist_mod == 0) {
			--r_num;
		}

		double time_1 = 0;
		double time_2 = 0;
		double wait_1 = 0;
		double wait_2 = 0;
		double delay_1 = 0;
		double delay_2 = 0;
		//求0 ...... cust_pos - c_iter......0的行驶时间time_1,等待时间wait_1,延迟时间delay_1
		time_1 = l_time[org_route.customers[cust_pos]];
		wait_1 = total_wait[org_route.customers[cust_pos]];
		delay_1 = total_delay[org_route.customers[cust_pos]];

		time_1 += peer_time[org_route.customers[cust_pos]][*c_iter];

		if (time_1 < customer[*c_iter].b_time) {					//小于最早时间
			wait_1 += customer[*c_iter].b_time - time_1;
			time_1 = customer[*c_iter].b_time;
		}
		else if (time_1 > customer[*c_iter].e_time)				//大于最迟时间
			delay_1 += time_1 - customer[*c_iter].e_time;
		//在时间窗之内
		time_1 += customer[*c_iter].s_time;
		for (vector<int>::iterator c_ind = c_iter + 1; c_ind != temp_route.customers.end(); ++c_ind) {
			time_1 += peer_time[*(c_ind - 1)][*c_ind];
			if (time_1 < customer[*c_ind].b_time) {
				wait_1 += customer[*c_ind].b_time - time_1;
				time_1 = customer[*c_ind].b_time;
			}
			else if (time_1>customer[*c_ind].e_time)
				delay_1 += time_1 - customer[*c_ind].e_time;
			time_1 += customer[*c_ind].s_time;
		}
		//求0 ...... c_iter - 1 -> cust_pos + 1......0的行驶时间time_2,等待时间wait_3,延迟时间delay_2
		time_2 = l_time[*(c_iter - 1)];
		wait_2 = total_wait[*(c_iter - 1)];
		delay_2 = total_delay[*(c_iter - 1)];

		time_2 += peer_time[*(c_iter - 1)][org_route.customers[cust_pos + 1]];
		if (time_2<customer[org_route.customers[cust_pos + 1]].b_time) {
			wait_2 += customer[org_route.customers[cust_pos + 1]].b_time - time_2;
			time_2 = customer[org_route.customers[cust_pos + 1]].b_time;
		}
		else if (time_2>customer[org_route.customers[cust_pos + 1]].e_time) delay_2 += time_2 - customer[org_route.customers[cust_pos + 1]].e_time;
		time_2 += customer[org_route.customers[cust_pos + 1]].s_time;
		for (int i = cust_pos + 2; i<org_route.customers.size(); ++i) {
			time_2 += peer_time[org_route.customers[i - 1]][org_route.customers[i]];
			if (time_2<customer[org_route.customers[i]].b_time) {
				wait_2 = customer[org_route.customers[i]].b_time - time_2;
				time_2 = customer[org_route.customers[i]].b_time;
			}
			else if (time_2>customer[org_route.customers[i]].e_time) delay_2 += customer[org_route.customers[i]].e_time - time_2;
			time_2 += customer[org_route.customers[i]].s_time;
		}

		if (r_lamda[0] == 1) {													//如果是优化f[0],用delta 记录 路径数	
			delta = r_num;
		}
		else if (r_lamda[1] == 1) {												//如果对f[1](总行驶距离)优化：
			delta = now_dist + dist_org + dist_mod;
		}
		else if (r_lamda[2] == 1) {                                            //如果对f[2]最长路径行驶时间优化（the longer time between both routes）
			if (time_2 > time_1)
				delta = time_2;
			else
				delta = time_1;

		}
		else if (r_lamda[3] == 1) {													//对f[3]总等待时间进行优化
			delta = now_wait + wait_1 + wait_2;
		}
		else {																	//对f[4]延迟时间优化
			delta = now_delay + delay_1 + delay_2;
		}

		if (min_delta > delta) {						//用min_dalta记录LS6领域操作后对应目标的值，其中对于f[2]，min_delta记录的是r_ind和i路径中最长路径行驶时间
			min_delta = delta;
			c_pos = c_iter - temp_route.customers.begin();						//c_pos是最优的插入位置
			feasible = true;													//返回为true，找到最优的位置了
		}

	}
	return feasible;

}

/*功能：将r_ind路径与解中其他路径进行LS6操作，
详情如下：
1.先找到最优的r_pos和min_c_pos
2.判断 是否找到最优插入位置或对f[2]的优化是否有意义
false	: 	将chrome.routes[r_ind]根据r_ind拆分成两个,如：
chrome.routes[temp_route_1] : 0,......,c_ind,0；
chrome.routes[temp_route_2] : 0,c_ind + 1,......,0
true	:	实现交换chrome.routes[r_ind] 和chrome.routes[r_pos]
交换前：
chrome.routes[r_ind] : 0,......,c_ind,c_ind + 1,......,0
chrome.routes[r_pos] : 0,......,min_c_pos,min_c_pos + 1,......,0
交换后：
chrome.routes[r_ind] : 0,......,c_ind,c_pos,min_c_pos,......,0
chrome.routes[r_pos] : 0,......,min_c_pos - 1,c_ind + 1,......,0
记录下最小的对应目标值。如：对于f[1],用delta记录最短的总行驶距离,用min_c_pos记录下最优的插入位置
r_pos：记录下进行与r_ind路径成功进行LS6操作的路径索引
r_ind:根据不同的目标选用不同策略先的路径
*/
void two_opt(Chromosome &chrome, double r_lamda[], int r_ind, double now_time) {

	//1.先找到最优的r_pos和min_c_pos
	int c_ind = rand() % (chrome.routes[r_ind].customers.size() - 1);
	compute_route_time(chrome.routes[r_ind]);
	compute_route_distance(chrome.routes[r_ind]);
	compute_route_capacity(chrome.routes[r_ind]);

	/*
	将r_ind路径与解中其他路径进行LS6操作，记录下最小的对应函数值。找到feasible = true;否则为false
	如：对于f[1],用delta记录最短的总行驶距离,用min_c_pos记录下最优的插入位置

	*/
	double dist = 0, now_wait, now_delay;
	int c_pos = 0;
	//记录使用LS6最优的插入位置索引
	int min_c_pos = 0;
	double min_delta = 99999999, delta = 0;
	bool feasible = false;
	int r_pos;									//记录下进行与r_ind路径成功进行LS6操作的路径索引
	for (int i = 0; i < chrome.routes.size(); ++i) {
		//求除r_ind和i 路径的路径长度之和、等待时间之和、延迟时间之和。
		if (i != r_ind) {
			//求除r_ind和i 路径的路径长度之和、等待时间之和、延迟时间之和。
			dist = 0; now_time = 0; now_wait = 0; now_delay = 0;
			for (int j = 0; j < chrome.routes.size(); ++j)
				if (j != r_ind && j != i) {
					dist += chrome.routes[j].travel_dist;				//除r_ind和i 路径的路径长度之和
					now_wait += chrome.routes[j].wait_time;				//除r_ind和i 路径的等待时间之和
					now_delay += chrome.routes[j].delay_time;			//除r_ind和i 路径的延迟时间之和
				}
			//用LS6,查找最优的插入位置：找到最优位置后，用delta记录对应目标的值，其中代表f[2]的表示：r_ind路径与i路径的最长行驶时间
			if (check_best_pos(chrome.routes[r_ind], c_ind, chrome.routes[i], c_pos, delta, r_lamda, chrome.routes.size(), dist, now_time, now_wait, now_delay)) {    //找到最优位置
				feasible = true;
				if (min_delta > delta) {
					min_delta = delta;
					min_c_pos = c_pos;
					r_pos = i;
				}
			}
		}
	}
	/*
	2.判断 是否找到最优插入位置或对f[2]的优化是否有意义
	false	: 	将chrome.routes[r_ind]根据r_ind拆分成两个,如：
	chrome.routes[temp_route_1] : 0,......,c_ind,0；
	chrome.routes[temp_route_2] : 0,c_ind + 1,......,0
	true	:	实现交换chrome.routes[r_ind] 和chrome.routes[r_pos]

	*/

	if (!feasible || (r_lamda[2] == 1 && min_delta > now_time)) {			//没找到最佳插入 || (对最长路径行驶时间f[2]进行优化 || LS6处理的f[2] > 原来的f[2] )
																			/*
																			将chrome.routes[r_ind]根据r_ind拆分成两个,如：
																			chrome.routes[temp_route_1] : 0,......,c_ind,0；
																			chrome.routes[temp_route_2] : 0,c_ind + 1,......,0
																			*/
		Route temp_route_1, temp_route_2;
		for (int i = 0; i <= c_ind; ++i)
			temp_route_1.customers.push_back(chrome.routes[r_ind].customers[i]);
		temp_route_1.customers.push_back(0);

		temp_route_2.customers.push_back(0);
		for (int i = c_ind + 1; i<chrome.routes[r_ind].customers.size(); ++i)
			temp_route_2.customers.push_back(chrome.routes[r_ind].customers[i]);
		compute_route_time(temp_route_1);
		compute_route_distance(temp_route_1);
		compute_route_time(temp_route_2);
		compute_route_distance(temp_route_2);

		//将路径chrome.routes[r_ind]从解中删除
		chrome.routes.erase(chrome.routes.begin() + r_ind);
		//将路径chrome.routes[temp_route_1]、chrome.routes[temp_route_2]中没有意义的路径也从解中删除
		if (temp_route_1.customers.size() != 2)  chrome.routes.push_back(temp_route_1);
		if (temp_route_2.customers.size() != 2)  chrome.routes.push_back(temp_route_2);
	}//找到了 ||  不是对f[2]进行优化 || (LS6处理的f[2] < 原来的f[2] )
	else {
		vector<int> r1, r2;
		//将0 ......c_ind 添加到r1中
		for (int i = 0; i <= c_ind; ++i)
			r1.push_back(chrome.routes[r_ind].customers[i]);
		//将(最优插入位置)min_c_pos ...... chrome.routes[r_pos].customers.size() -1 添加到r1中
		for (int i = min_c_pos; i<chrome.routes[r_pos].customers.size(); ++i)
			r1.push_back(chrome.routes[r_pos].customers[i]);
		//将0 ......min_c_pos-1 添加到r2中
		for (int i = 0; i<min_c_pos; ++i)
			r2.push_back(chrome.routes[r_pos].customers[i]);
		//将c_ind+1 ......chrome.routes[r_ind].customers.size() 添加到r2中
		for (int i = c_ind + 1; i<chrome.routes[r_ind].customers.size(); ++i)
			r2.push_back(chrome.routes[r_ind].customers[i]);
		//分别将客户点容器赋给对应路径的客户点容器变量
		chrome.routes[r_ind].customers = r1;
		chrome.routes[r_pos].customers = r2;
		//分别计算索引为r_ind和r_pos路径的行驶时间、行驶距离等；
		compute_route_time(chrome.routes[r_ind]);
		compute_route_distance(chrome.routes[r_ind]);
		compute_route_time(chrome.routes[r_pos]);
		compute_route_distance(chrome.routes[r_pos]);
		//如果路径变为无效路径，从解中删除此路径
		if (r1.size() == 2) chrome.routes.erase(chrome.routes.begin() + r_ind);
		if (r2.size() == 2) chrome.routes.erase(chrome.routes.begin() + r_pos);
	}
}

/*功能：		检测将cust插入r_iter路径中是否是最优
问题：（706，717）（delay_time += customer[*c_iter].e_time - travel_time;）原来是这样的，这样可以算出延迟时间吗？
check_best_pos(cust,*r_iter,c_pos,delta,r_lamda,now_dist,now_time,now_wait,now_delay)
cust_index	:	cust			客户点在总客户点容器中的索引
route		:	r_iter			解中要插入索引为cust客户点的路径
c_pos		:	c_pos			记录cust_indx插入到route中最佳的位置索引
min_delta	:	delta			插入后，对应目标的函数值，f[2]除外；
now_dist	:	now_dist		没有插入custermer[cust]客户点时，解的总行驶距离
now_time	:	now_time		没有插入custermer[cust]客户点时，解的总行驶时间
now_wait	:	now_wait		除r_iter以外其他路径的等待时间之和
now_delay	:	now_delay		除r_iter以外其他路径的延迟时间之和
*/
bool check_best_pos(int cust_index, Route& route, int &c_pos, double &min_delta, double r_lamda[], double now_dist, double now_time, double now_wait, double now_delay) {

	bool feasible = false;
	min_delta = 99999999; double delta_dist = 0, delta_time = 0;
	double capacity = 0, delta = 0;
	//容量约束检查
	vector<int>::iterator c_iter;
	for (c_iter = route.customers.begin() + 1; c_iter != route.customers.end() - 1; ++c_iter) {
		capacity += customer[*c_iter].demand;
	}
	capacity += customer[cust_index].demand;
	if (capacity>max_capacity)					//容量约束没通过
		return false;
	//容量约束检查通过，检查时间约束
	compute_route_time(route);
	for (c_iter = route.customers.begin() + 1; c_iter != route.customers.end(); ++c_iter) {
		double travel_time = 0;
		double wait_time = total_wait[*(c_iter - 1)];
		double delay_time = total_delay[*(c_iter - 1)];
		//插入到c_iter之前，求到达cust_index的时间
		a_time[cust_index] = l_time[*(c_iter - 1)] + peer_time[*(c_iter - 1)][cust_index];
		//如果到达cust_index的时间未通过时间约束，继续for循环
		if (a_time[cust_index] > customer[cust_index].e_time + max_delay_time)
			continue;
		//提前到cust_index客户点
		if (a_time[cust_index] < customer[cust_index].b_time) {
			wait_time += customer[cust_index].b_time - a_time[cust_index];
			a_time[cust_index] = customer[cust_index].b_time;
		}
		//迟到cust_index客户点，但迟到时间不超过max_delay_time
		if (a_time[cust_index] > customer[cust_index].e_time)
			delay_time += a_time[cust_index] - customer[cust_index].e_time;

		l_time[cust_index] = a_time[cust_index] + customer[cust_index].s_time;
		//求插入cust_index的c_iterr 增量时间,即：插入cust_index客户点后，车辆到达客户点i的时间减去插入前车辆到过客户点i的时间
		delta_time = l_time[cust_index] + peer_time[cust_index][*c_iter] - a_time[*c_iter];
		//插入cust_index的c_iterr 增量时间 > c_iter的最大合法延迟时间
		if (delta_time>max_wait[*c_iter])
			continue;
		//在c_iter的距离增量
		delta_dist = -peer_distance[*(c_iter - 1)][*c_iter] + peer_distance[*(c_iter - 1)][cust_index] + peer_distance[cust_index][*c_iter];
		//到c_iter的行驶时间
		travel_time = l_time[cust_index] + peer_time[cust_index][*c_iter];
		//早到c_iter
		if (travel_time < customer[*c_iter].b_time) {
			wait_time += customer[*c_iter].b_time - travel_time;
			travel_time = customer[*c_iter].b_time;
		}
		//迟到c_iter
		if (travel_time > customer[*c_iter].e_time)
			//delay_time += customer[*c_iter].e_time - travel_time;						//原来是这样的，这样可以算出延迟时间吗
			delay_time += travel_time - customer[*c_iter].e_time;

		travel_time += customer[*c_iter].s_time;
		for (vector<int>::iterator c_ind = c_iter; c_ind != route.customers.end() - 1; ++c_ind) {
			travel_time += peer_time[*c_iter][*(c_iter + 1)];
			if (travel_time < customer[*(c_iter + 1)].b_time) {
				wait_time += customer[*(c_iter + 1)].b_time - travel_time;
				travel_time = customer[*(c_iter + 1)].b_time;
			}
			//travel_time表示到达c_iter + 1的行驶时间，
			if (travel_time > customer[*(c_iter + 1)].e_time)
				//delay_time+=customer[*(c_iter+1)].e_time-travel_time;  //老师给的代码是这样的
				delay_time += travel_time - customer[*(c_iter + 1)].e_time; //我认为是这样的

			travel_time += customer[*(c_iter + 1)].s_time;
		}

		if (r_lamda[1] == 1) {
			delta = now_dist + delta_dist;
		}
		else if (r_lamda[2] == 1) {
			delta = travel_time;
		}
		else if (r_lamda[3] == 1) {
			delta = now_wait + wait_time;
		}
		else
			delta = now_delay + delay_time;

		if (min_delta>delta) {
			min_delta = delta;
			c_pos = c_iter - route.customers.begin();
			feasible = true;
		}

	}
	return feasible;

}

/*功能:判断cust客户点是否可以插入到route路径中,并且实现插入
1.route			:	某路径
2.cust			:	路径中的客户点
*/
bool check_feasible_place(Route &route, int cust) {
	double time = 0, capacity = 0;
	vector<int>::iterator c_iter;
	compute_route_time(route);
	//检查添加客户点cust后的容量约束
	for (c_iter = route.customers.begin() + 1; c_iter != route.customers.end() - 1; ++c_iter) {
		capacity += customer[*c_iter].demand;
	}
	capacity += customer[cust].demand;
	if (capacity > max_capacity)							//没通过容量约束
		return false;
	//时间约束检查：检查route中有没有合适插入的位置
	for (c_iter = route.customers.begin() + 1; c_iter != route.customers.end() - 1; ++c_iter) {					//将cust插入到c_iter所在的位置
		double tmp_time = time + peer_time[*(c_iter - 1)][cust];
		if (tmp_time > customer[cust].e_time + max_delay_time) {													//插入,c_iter之前的客户点有问题，不能插入,采取一般处理,继续往下搜
			time += peer_time[*(c_iter - 1)][*c_iter];
			if (time < customer[*c_iter].b_time) time = customer[*c_iter].b_time;
			time += customer[*c_iter].s_time;
			continue;
		}
		//可以插入,c_iter之前的客户点没问题
		if (tmp_time<customer[cust].b_time) tmp_time = customer[cust].b_time;										//将tmp_time设置为cust客户点的b_time
																													//判断c_iter之后的客户点怎样
		tmp_time += customer[cust].s_time + peer_time[cust][*c_iter];
		if (tmp_time <= a_time[*c_iter] || tmp_time <= customer[*c_iter].b_time) {								//判断c_iter后面的客户点也可以插入
																												//这里插入的是第一个最合理的位置
			route.customers.insert(c_iter, cust);																//将客户点cust插入到c_iter所在的位置,c_iter客户点后移.
			compute_route_capacity(route);
			compute_route_distance(route);
			compute_route_time(route);
			return true;																						//插入到合适的位置,退出函数
		}//如果c_iter后面的客户点不可以插入,继续正常操作
		time += peer_time[*(c_iter - 1)][*c_iter];
		if (time<customer[*c_iter].b_time) time = customer[*c_iter].b_time;
		time += customer[*c_iter].s_time;
	}//route路径中[0,route.customers.end() - 2]没有合适的位置,插入到0之前的位置
	time += peer_time[*(c_iter - 1)][cust];
	if (time > customer[cust].e_time + max_delay_time)																//不能插入
		return false;
	if (time<customer[cust].b_time)
		time = customer[cust].b_time;
	time += customer[cust].s_time + peer_time[cust][*c_iter];
	if (time <= customer[*c_iter].e_time) {																			//cust可以插入到0之前的位置
		route.customers.insert(c_iter, cust);
		compute_route_capacity(route);
		compute_route_distance(route);
		compute_route_time(route);
		return true;
	}//cust不可以插入到0之前的位置
	return false;

}

void check(Chromosome &chrome);

/*删除路径route_ind
将route_ind路径上每个的客户点尝试插入其他路径，最后只剩下0，0客户点时，删除路径
删除成功，返回true
其他返回false
*/
bool route_number(Chromosome &chrome, int route_ind) {
	vector<int>::iterator c_iter;
	vector<Route>::iterator r_iter;

	//对route_ind路径中的客户点进行遍历
	for (c_iter = chrome.routes[route_ind].customers.begin() + 1; c_iter != chrome.routes[route_ind].customers.end() - 1;) {
		//计算路径中各个客户点的行驶时间、等待时间、延迟时间和最大合法延迟时间（方便插入客户点）
		compute_route_time(chrome.routes[route_ind]);
		/*
		路径route_ind中可以删除的客户点要满足的条件:
		删除后:
		1.当删除点的下一客户点是0(车场):删除后,到达车场的时间要小于车场的最迟关闭时间
		2.当删除的下一个客户点不是车场:删除后,到达下一个客户点的时间小于下个客户点的最迟时间+最大允许迟到时间
		*/
		if ((*(c_iter + 1) == 0 && l_time[*(c_iter - 1)] + peer_time[*(c_iter - 1)][*(c_iter + 1)] < customer[0].e_time) ||								//删除的下一个客户点是车场
																																						//2.c_iter != 0	
			(*(c_iter + 1) != 0 && l_time[*(c_iter - 1)] + peer_time[*(c_iter - 1)][*(c_iter + 1)] < customer[*(c_iter + 1)].e_time + max_delay_time)) {		//删除的下一个客户点不是车场
			for (r_iter = chrome.routes.begin(); r_iter != chrome.routes.end(); ++r_iter) {													//路径0 - n-2,不包含n - 1(最后一条路径)													
				if (r_iter - chrome.routes.begin() != route_ind && check_feasible_place(*r_iter, *c_iter)) {								//*c_iter客户点可以插入到*r_iter路径
					c_iter = chrome.routes[route_ind].customers.erase(c_iter);//因为c_iter客户点被移除，c_iter = chrome.routes[route_ind].customers.begin()+1是容器中其它客户点
					break;
				}
			}
			if (r_iter == chrome.routes.end())
				return false;
		}
		else
			return false;

		if (r_iter == chrome.routes.end())																		//c_iter都没插入，到最后一条路径 ，不问它了，插入下一个
			++c_iter;																												//循环步长
	}
	if (chrome.routes[route_ind].customers.size() == 2) {
		chrome.routes.erase(chrome.routes.begin() + route_ind);																		//删除route_ind路径
		return true;
	}
	return false;
}

/*	此处有问题
向解chrome中插入单个客户点cust
insert_customer(chrome,r_lamda, customers[iter],now_time);

*/
void insert_customer(Chromosome &chrome, double r_lamda[], int cust, double now_time = 9999999) {
	//没有插入custermer[cust]客户点时，解的总行驶距离
	double now_dist = 0;
	//没有插入custermer[cust]客户点时，解的总等待时间
	double now_wait = 0;
	////没有插入custermer[cust]客户点时，解的总延迟时间
	double now_delay = 0;

	Route temp_route;
	vector<Route>::iterator r_iter, r_pos;
	//求目前总行驶距离
	for (r_iter = chrome.routes.begin(); r_iter != chrome.routes.end(); ++r_iter) {
		now_dist += r_iter->travel_dist;
	}

	int c_pos = 0;
	int min_c_pos = 0;
	double min_delta = 99999999;/*double min_time=0;*/
	double delta = 0;
	bool feasible = false;

	for (r_iter = chrome.routes.begin(); r_iter != chrome.routes.end(); ++r_iter) {
		now_wait = 0;
		now_delay = 0;
		for (vector<Route>::iterator r_index = chrome.routes.begin(); r_index != chrome.routes.end(); ++r_index)
			if (r_iter != r_index) {
				//除r_iter以外其他路径的等待时间
				now_wait += r_index->wait_time;
				//除r_iter以外其他路径的延迟时间
				now_delay += r_index->delay_time;
			}
		//此处有问题
		if (check_best_pos(cust, *r_iter, c_pos, delta, r_lamda, now_dist, now_time, now_wait, now_delay)) {  //检测将cust插入r_iter路径中是否是最优
			feasible = true;
			if (min_delta>delta) {
				min_delta = delta;
				min_c_pos = c_pos;
				r_pos = r_iter;
			}
		}
	}

	//没有找到最佳插入位置
	if (!feasible /*|| delta< min_delta*/ || (r_lamda[2] == 1 && now_time<min_delta)) {
		//新建一个路径，0,cust,0
		clear_route(temp_route);
		temp_route.customers.push_back(0);
		temp_route.customers.push_back(cust);
		temp_route.customers.push_back(0);

		compute_route_time(temp_route);
		compute_route_distance(temp_route);
		compute_route_capacity(temp_route);
		chrome.routes.push_back(temp_route);
	}
	else { //找到最佳插入位置，把cust插入到最佳路径r_pos的最佳位置min_c_pos中		 
		r_pos->customers.insert(r_pos->customers.begin() + min_c_pos, cust);
		compute_route_time(*r_pos);
		compute_route_distance(*r_pos);
	}
}

bool update_EP(Chromosome &new_chrome);
bool is_equal(Chromosome &chrome_1, Chromosome &chrome_2);

/*原来代码没有注释
void get_RP(int RP[],int size){

for(int i=0;i<size;++i)
RP[i]=i;
for(int i=0;i<size;++i){
int r0=i+rand()%(size-i);
std::swap(RP[i],RP[r0]);
}
}
*/
int RP[MAX];

/*
之前操作：从route中随机删除count个客户点 950
检查解中此路径是否有效：
有效返回：true
无效返回：false
//只检查这个路径上的时间约束，没有检查容量约束

问题：只检查时间约束，不检查容量约束，会不会删除的是取货点，车辆的货物不减少，会不会造成后面的车的货物超载
*/
bool check_feasible_route(Route &route) {
	double time = 0;
	//只检查了时间约束，没有检查容量约束
	for (vector<int>::iterator c_iter = route.customers.begin(); c_iter != route.customers.end() - 1; ++c_iter) {
		time += peer_time[*c_iter][*(c_iter + 1)];
		//不是有效路径的两种情况
		if (*(c_iter + 1) == 0 && time>customer[*(c_iter + 1)].e_time)			//如果下一个路径是0，若行驶时间 > 下一个客户点（车场）的最后时间，返回false.（即：不是有效路径）
			return false;
		//如果下一个路径不是0，若行驶时间 > 下一个客户点的关闭时间 + 最大延迟时间，返回false.（即：不是有效路径）
		if (*(c_iter + 1) != 0 && time>customer[*(c_iter + 1)].e_time + max_delay_time)
			return false;
		//到目前为止，是有效路径，且行驶时间 < 下一个客户点的最早开放时间
		if (time < customer[*(c_iter + 1)].b_time)
			time = customer[*(c_iter + 1)].b_time;
		time += customer[*(c_iter + 1)].s_time;
	}
	return true;
}

/*
功能：拆分解chrome的第r_ind条路径.产生新的路径添加到解的路径集合中
有问题:最终的temp_route中少一个连接一个0
*/
void split_route(Chromosome &chrome, vector<int> &tempCusts, int r_ind) {
	Route temp_route;
	//获得要拆分的路径
	vector<Route>::iterator r_iter = chrome.routes.begin() + r_ind;
	//
	vector<int>::iterator c_iter;
	double time = 0;
	temp_route.customers.push_back(0);
	//这是索引还是客户点
	int previous_cust = 0, now_cust, tmp = 0, tmp2 = 0;

	for (c_iter = r_iter->customers.begin(); c_iter != r_iter->customers.end() - 1; ++c_iter) {
		now_cust = *(c_iter + 1);
		time += peer_time[previous_cust][now_cust];
		double return_time = time;

		if (return_time<customer[now_cust].b_time)
			return_time = customer[now_cust].b_time;
		return_time += customer[now_cust].s_time + peer_time[now_cust][0];

		if (time > customer[now_cust].e_time + max_delay_time || return_time>customer[0].e_time) {
			temp_route.customers.push_back(0);
			if (temp_route.customers.size() != 2)
			{
				compute_route_time(temp_route);
				compute_route_distance(temp_route);
				compute_route_capacity(temp_route);
				chrome.routes.push_back(temp_route);
			}
			temp_route.customers.clear();

			temp_route.customers.push_back(0);
			previous_cust = 0;
			time = 0;
			tmp = *c_iter;
			//如果第1个客户点就不满足条件,怎么办,如 0,6
			if (tmp <= 0 && now_cust > 0)
			{
				tempCusts.push_back(now_cust);
			}
		}
		else {
			temp_route.customers.push_back(now_cust);
			if (time<customer[now_cust].b_time) time = customer[now_cust].b_time;
			time += customer[now_cust].s_time;
			previous_cust = now_cust;
		}
	}
	if (temp_route.customers.size() != 2)
	{
		compute_route_time(temp_route);
		compute_route_distance(temp_route);
		compute_route_capacity(temp_route);
		chrome.routes.push_back(temp_route);
	}
	chrome.routes.erase(chrome.routes.begin() + r_ind);
}

int adaptive_select_neighbor(double r_lamda[])
{
	int ind = 0;							//记录要优化的目标
	vector <double> ls;
	for (int i = 0; i < FUNC_NUM; i++)
		if (r_lamda[i] == 1)
		{
			ind = i;
			break;
		}
	//假如采用轮盘赌的方法，选择目标索引
	double sum = 0;
	for (int i = 0; i < NEIGH_NUM; i++)
	{
		ls.push_back(LSArray[ind][i]);
		sum += LSArray[ind][i];
	}
	for (int i = 0; i < NEIGH_NUM; i++)
		ls[i] /= sum;
	ind2 = RouletteWheelSelection(ls);		//选择第ind2个局部搜索
	LSTable[ind2] = true;
	lsn[ind2]++;
	return ind2;
}

//更新目标-局部搜索矩阵
void updateobjls()
{
	double s[FUNC_NUM][NEIGH_NUM];
	//每个目标提升度的总和
	double objsum[FUNC_NUM];
	//计算每个局部搜索对每个目标的平均提升度
	for (int i = 0; i < FUNC_NUM; i++)						//第i个目标在各个局部搜索中的提升
	{
		objsum[i] = 0;
		for (int j = 0; j < NEIGH_NUM; j++)					//第j个领域搜索
		{
			s[i][j] = FIR[i][j] / lsn[j];
			objsum[i] += s[i][j];
			//将提升度归为0
			FIR[i][j] = 0;
		}
	}//规一化数组

	for (int i = 0; i < FUNC_NUM; i++)						//第i个目标在各个局部搜索中的提升
		for (int j = 0; j < NEIGH_NUM; j++)					//第i个目标在第j个局部搜索上的提升度
			if (objsum[i] != 0)
				s[i][j] /= objsum[i];						//相当于CVi，求提升的平均值
			else
				s[i][j] = 0;

	//质量更新
	for (int i = 0; i < FUNC_NUM; i++)
		for (int j = 0; j < NEIGH_NUM; j++)
			q[i][j] = 0.5 * q[i][j] + 0.5 * s[i][j];
	//可能性更新
	for (int i = 0; i < FUNC_NUM; i++)			//目标的遍历
	{
		double pmin = LSArray[i][0], qsum = 0;
		for (int j = 0; j < NEIGH_NUM; j++)    //找到此目标对应局部搜索中，最小的。
		{
			if (pmin > LSArray[i][j])
				pmin = LSArray[i][j];
			qsum += q[i][j];
		}

		objsum[i] = 0;
		for (int j = 0; j < NEIGH_NUM; j++)
		{
			if (qsum != 0)
				LSArray[i][j] = pmin + (1 - NEIGH_NUM * pmin) * (q[i][j] / qsum);
			else
				LSArray[i][j] = pmin;
			objsum[i] += LSArray[i][j];
		}

		//目标对应局部搜索规一化
		for (int j = 0; j < NEIGH_NUM; j++)
			LSArray[i][j] /= objsum[i];
	}
	//标记
	for (int i = 0; i < NEIGH_NUM; i++)
	{
		LSTable[i] = false;
		lsn[i] = 0;
	}
	is_EP = false;
}

/*
功能：remove k customers and re-inserted them to the best place one by one
采用三种不同的局部搜索方法{"删除单个客户点","删除n个客户点","用LS6领域操作算子"}，对解进行优化，将新解进行执行update_EP(chrome);操作
有问题:最终的temp_route中少一个连接一个0 ,			split_route(chrome,r_ind);
有问题（821）：只检查时间约束，不检查容量约束，会不会删除的是取货点，车辆的货物不减少，会不会造成后面的车的货物超载		check_feasible_route(chrome.routes[r_ind])
*/
void localsearch_k(Chromosome &chrome, int &count, double r_lamda[], Chromosome &best, int index) {

	//如果LS都选择一遍后，且有解加入到EP中，则更新table
	if (LSTable[0] && LSTable[1] && LSTable[2] && is_EP)
	{
		//每个局部搜索在各个目标上提升度的平均值
		updateobjls();
	}

	Chromosome previous = chrome;
	//存储的是客户点
	vector<int> customers;					//删除客户点的集合,存储的是客户点在全部客户点容器中的索引
	Route temp_route;
	double now_time = 0;						//记录一条路径的行驶时间,当优化f[1]、f[3]、f[4]，now_time = 0;当优化f[2]，now_time = 总行驶时间最长的路径的行驶时间
	customers.clear();
	/*功能：不同目标上的选择路径操作
	从解中选择路径，其中最长路径行驶时间f[3]：它选择路径的策略：选择总行驶时间最长的路径；
	f[2]、f[4]、f[5]选择路径的策略：用轮盘赌选择策略，选择较长路径
	*/
	int r_ind = 0;

	if (r_lamda[2] != 1)															///f[2]、f[4]、f[5]选择路径的策略
		r_ind = rand() % chrome.routes.size();												//随机选择一条路径,now_time = 0;
	else {																		//f[3]（最长路径行驶时间）：选择路径的策略，选择总行驶时间最长的路径；
		for (int i = 0; i<chrome.routes.size(); ++i)											//选择总行驶时间最长的路径 
			if (now_time<chrome.routes[i].travel_time) {
				now_time = chrome.routes[i].travel_time;
				r_ind = i;
			}
	}
	/*
	{"删除单个客户点","删除n个客户点","用LS6领域操作算子"}，对解进行改造。选中的概率均为1/3.
	*/
	//此处的随机选择，改成自适应的选择
	//double p=rand()%1000/1000.0;				//生成[0,1]的随机数，这样可以生成精确度到千分位上的[0,1]小数
	int k = adaptive_select_neighbor(r_lamda);
	//if(p < 1.0/3 ){								//删除1个随机的客户点
	if (k == 0) {
		//		int r_ind=r_ind_1;
		int cust_ind = rand() % (chrome.routes[r_ind].customers.size() - 2) + 1;		//随机的求客户点的索引
		int cust = chrome.routes[r_ind].customers[cust_ind];						//得到r_ind路径上的索引为cust_ind的客户点
		customers.push_back(cust);												//将客户点插入删除集合
		chrome.routes[r_ind].customers.erase(chrome.routes[r_ind].customers.begin() + cust_ind);		//从路径中删除客户点
		if (chrome.routes[r_ind].customers.size() == 2)							//如果路径中没有有用客户点，则从解中删除此路径
			chrome.routes.erase(chrome.routes.begin() + r_ind);
		else {//是有意义的路径
			  //删除chrome.routes[r_ind]路径的cust_ind后，检查删除后的路径是否有效	
			if (check_feasible_route(chrome.routes[r_ind])) {//有问题：只检查时间约束，不检查容量约束，会不会删除的是取货点，车辆的货物不减少，会不会造成后面的车的货物超载
				compute_route_time(chrome.routes[r_ind]);
				compute_route_distance(chrome.routes[r_ind]);
			}
			else {																//路线无效,就拆分此路线
																				//有问题的路径:最终的temp_route中少一个连接一个0
				split_route(chrome, customers, r_ind);										//拆分路径,产生新的路径加入到解中
			}
		}
	}
	//else if(p < 2.0/3 && p >= 1.0/3 ) {		//删除随机数count个随机的客户点	
	else if (k == 1) {
		int count = rand() % (chrome.routes[r_ind].customers.size() - 2) + 1;			//随机要删除的客户点数在[1,chrome.routes[r_ind].customers.size()-2]之间
		for (int i = 0; i < count; ++i) {
			int cust_ind = rand() % (chrome.routes[r_ind].customers.size() - 2) + 1;	//随机产生客户点索引
			customers.push_back(chrome.routes[r_ind].customers[cust_ind]);		//将该客户点插入到删除容器
			chrome.routes[r_ind].customers.erase(chrome.routes[r_ind].customers.begin() + cust_ind);	//从解中的路径中删除此客户点
		}
		if (chrome.routes[r_ind].customers.size() == 2)							//如果解中的这条路径无效,从解中删该路径
			chrome.routes.erase(chrome.routes.begin() + r_ind);
		else {		//路径长  > 2		问题：只检查时间约束，不检查容量约束，会不会删除的是取货点，车辆的货物不减少，会不会造成后面的车的货物超载
			if (check_feasible_route(chrome.routes[r_ind])) {
				compute_route_time(chrome.routes[r_ind]);
				compute_route_distance(chrome.routes[r_ind]);
			}
			else {
				split_route(chrome, customers, r_ind);
			}
		}
	}
	else	//使用LS6领域操作产生新解，不删除客户点，由于传入的r_lamda[]数组不同，优化的方向也不同
		two_opt(chrome, r_lamda, r_ind, now_time);


	//将删除容器中的客户点插入到解中产生新解
	for (int iter = 0; iter < customers.size(); ++iter)
		insert_customer(chrome, r_lamda, customers[iter], now_time);		//向解chrome中插入单个客户点customers[iter]
	compute_f(chrome);								//chrome在局部搜索中变化
	compute_f(best);								//best记录下原来的chrome的值
	int i;
	for (i = 0; i<FUNC_NUM; ++i)
		if (r_lamda[i] == 1)
			break;
	if (i != FUNC_NUM) index = i;
	if (best.f[index] > chrome.f[index]) {
		best = chrome;
		//		count=0;
	}
	else {
		if (chrome.f[index] > best.f[index]) chrome = best;
		//	++count;
	}
	/*}*/
	++count;
	if (update_EP(chrome))
	{
		double div = 0;
		for (int i = 0; i < FUNC_NUM; i++)       //目标
		{
			//df[i] = 0;
			div = previous.f[i] - chrome.f[i];
			if (div > 0)
				FIR[i][ind2] += (1.0 * div) / previous.f[i];					//局部搜索ind2对目标i的提升累加
																				//sum += df[i];
		}
		is_EP = true;
	}
	//if(!update_EP(chrome))
	//BETWEEN.push_back(chrome);
}

/*功能 ：对解进行局部搜索，产生新解，更新存档，对不同的目标，局部搜索不同
分别采用
删除单个路径   ，优化f[0];
或
随机选取{"删除单个客户点","删除n个客户点","用LS6领域操作算子"}，优化f[1]、f[2]、f[3]、f[4]
最终都要更新存档EP
*/
void local_search(Chromosome &chrome, double r_lamda[], int index) {
	if (r_lamda[0] == 1) {			//通过删除单个路径，来局部搜索，优化第一个目标
									/*
									选择解中客户点数最少的路径，删除，将其客户点分别插入解的其它路径
									*/
		while (true) {
			//定义一个迭代器类型
			vector<int>::iterator c_iter;
			vector<Route>::iterator r_iter;
			int route_ind = 0;																			//记录第几条路径
																										//找到解中找客户点最少的路径（即：简单的轮盘赌选择策略），记录其下标，尝试删除路径。
			int min_num = 999999;
			for (r_iter = chrome.routes.begin(); r_iter != chrome.routes.end(); ++r_iter) {
				if (min_num > r_iter->customers.size()) {
					min_num = r_iter->customers.size();
					route_ind = r_iter - chrome.routes.begin();
				}
			}
			if (!route_number(chrome, route_ind))														//删除route_ind路径失败，跳出while循环
				break;
			//如果删除路径成功，继续进行删除
		}
		//对新解的所有路径，求其行驶路径、容量、行驶时间
		for (vector<Route>::iterator iter = chrome.routes.begin(); iter != chrome.routes.end(); ++iter) {
			compute_route_distance(*iter);
			compute_route_capacity(*iter);
			compute_route_time(*iter);
		}
		//计算新解的是目标函数值
		compute_f(chrome);																		//更新新解的5个目标函数值
																								//把新解插入存档和边界解，并优化存档	
		update_EP(chrome);
		//if(!update_EP(chrome))
		//BETWEEN.push_back(chrome);
	}
	else {																						//分别对f[1]、f[2]、f[3]、f[4]进行优化																					//
		int count = 0; Chromosome best = chrome;
		for (count = 0; count<1;) {
			//采用三种不同的局部搜索方法{"删除单个客户点","删除n个客户点","用LS6领域操作算子"}，对解进行优化，将新解进行执行cout++;update_EP(chrome);操作
			localsearch_k(chrome, count, r_lamda, best, -1);
		}
		for (vector<Route>::iterator iter = chrome.routes.begin(); iter != chrome.routes.end(); ++iter) {
			compute_route_distance(*iter);
			compute_route_capacity(*iter);
			compute_route_time(*iter);
		}
		compute_f(chrome);
		//	update_EP(chrome,-1);
	}

}

/*判断存档中的解与新解是否相等
chrome_1:存档中的解；
chrome_2:新解
只要两个解的某一方向的值大于0.00001，两个解就不相等。返回false
其它返回true
*/
bool is_equal(Chromosome &chrome_1, Chromosome &chrome_2) {

	for (int i = 0; i<FUNC_NUM; ++i) {
		if (fabs(chrome_1.f[i] - chrome_2.f[i])>0.000001) return false;			//只要两个解的某一方向的值大于0.00001，两个解就不相等。返回false
	}
	return true;
}

/*
如果chrome_1与chrome2的关联矩阵都相等，返回true
否则，返回false
*/
bool is_equal_box(Chromosome &chrome_1, Chromosome &chrome_2) {

	for (int i = 0; i<FUNC_NUM; ++i) {
		if (abs(chrome_1.box_f[i] - chrome_2.box_f[i])>0.000001) return false;
	}
	return true;
}
//计算存档A中每个非占优解对应的一个关联矩阵B = {B0,B1,B2,B3,B4}
void get_box(Chromosome &chrome) {
	for (int i = 0; i<FUNC_NUM; ++i) {
		if (chrome.f[i] == 0) chrome.box_f[i] = 0;
		else chrome.box_f[i] = log(chrome.f[i]) / log(1 + epsilon);					//此处有问题，论文中是：chrome.box_f[i]=log(chrome.f[i] + 1)/log(1+epsilon);
	}
}

//在某个方向上好，就返回true
bool is_extreme(Chromosome &chrome) {
	for (int index = 0; index<FUNC_NUM; ++index) {
		if (extreme[index].f[index] - chrome.f[index]>0.00001) return true;//在这个方向上好，就返回true
	}
	return false;
}

//存储边界解：从5个目标考查新解，将新解存储在最优目标上，值越小越好
void preserve_extreme(Chromosome &chrome, int index) {

	if (extreme[index].f[index] - chrome.f[index]>0.00001) {				//在这个方向上好，得到保存，极限
		extreme[index] = chrome;
	}

}


//更新边界解
void update_extreme(Chromosome &new_chromosome) {

	for (int i = 0; i<FUNC_NUM; ++i) {					//对5个目标进行考查
		preserve_extreme(new_chromosome, i);
	}
}

/*
功能：把新解插入存档和边界解（可能插入，也可能不插入），并优化存档
更新新解new_chromosome
1.如果存档为空
更新边界解
向存档中添加新解
2.存档不为空
if new chromosome is worse than the one in EP, return，不更新存档

否则，更新存档
判断新解是否是极值解
是边界解
将新解与存档中的解比较，如果存档中的某个解没有优于新解，就删除存档中的那个解
最后将新解插入存档
不是边界解
将新解与存档中的解比较
如果新解优于存档解，将新解加入存档
如果存档中的某个解没有优于新解，就删除存档中的那个解
新解比所有存档中的解都差，且不是极值解
将新解的关联矩阵与存档解的关联矩阵进行比较
存档中的一个解的关联矩阵与新解的相等，退出循环，不再比较
存档中的解的关联矩阵与新解的都不相等，且存档的解也不优于新解。		将新解插入存档

*/
bool update_EP(Chromosome &new_chromosome) {

	bool flag = false, f = false;
	//计算存档A中每个解的关联矩阵
	get_box(new_chromosome);

	if (EP.size() == 0) {								//如果解集为空
														//更新边界解
		update_extreme(new_chromosome);
		//向存档中添加新解
		EP.push_back(new_chromosome);				//把新解加入存档中
		f = true;
	}
	else {											//解集不为空
		vector<Chromosome>::iterator iter = EP.begin();
		for (; iter != EP.end(); ++iter) {                            // if new chromosome is worse than the one in EP, return
			if (is_better(*iter, new_chromosome) || is_equal(*iter, new_chromosome)) return false;
		}
		iter = EP.begin();
		//判断新解是否是极值解
		if (is_extreme(new_chromosome)) {					//	如果新解是极值解
			update_extreme(new_chromosome);
			iter = EP.begin();
			while (iter != EP.end()) {
				if (is_better(new_chromosome, *iter)) //将新解与存档中的解比较，如果存档中的某个解没有优于新解，就删除存档中的那个解
					iter = EP.erase(iter);
				else ++iter;
			}
			EP.push_back(new_chromosome);				//最后将新解插入存档
			f = true;

		}
		else {											//新解不是极值解
			iter = EP.begin();
			while (iter != EP.end()) {
				if (is_better(new_chromosome, *iter)) { //将新解与存档中的解比较，如果存档中的某个解没有优于新解，就删除存档中的那个解
					iter = EP.erase(iter);
					flag = true;
				}
				else ++iter;
			}
			if (flag) {									 //如果新解优于存档解，将新解加入存档
				EP.push_back(new_chromosome);
				f = true;
			}
			else {										//新解比所有存档中的解都差，且不是极值解
				iter = EP.begin();
				for (; iter != EP.end(); ++iter) {			//将新解的关联矩阵与存档解的关联矩阵进行比较
					if (is_equal_box(new_chromosome, *iter) || is_better_box(*iter, new_chromosome)) break;
				}										//存档中的解的关联矩阵与新解的都不相等，且存档的解也不优于新解
				if (iter == EP.end()) {
					EP.push_back(new_chromosome);		//将新解插入存档
					f = true;
				}
			}
		}
	}
	return f;
}


/*
功能：判断new_chromosome是否是占优解，删除其中相对new_chromosome不占优的，将占优的new_chromosome解加入占优解集total_best中
若new_chromosome加入占优解集total_best成功，返回true
否则，加入未成功，返回false
*/

bool update_best(Chromosome &new_chromosome) {

	int size = EP.size();
	if (total_best.size() == 0)
		total_best.push_back(new_chromosome);
	else {
		//将存档EP中的解与new_chromosome进行比较
		vector<Chromosome>::iterator iter = total_best.begin();
		for (; iter != total_best.end(); ++iter) {
			//将容器中的所有iter与new_chromosome比较，iter是否优于new_chromosome    ||     iter与chomosome是否相等
			if (is_better(*iter, new_chromosome) || is_equal(*iter, new_chromosome))
				return false;
		}//new_chromosome优于EP中所有的解

		 //从EP中去除劣于new_chromosome的解
		iter = total_best.begin();
		while (iter != total_best.end()) {
			if (is_better(new_chromosome, *iter)) iter = total_best.erase(iter);
			else ++iter;
		}
		//将new_chromosome加入占优解集
		total_best.push_back(new_chromosome);
	}
	return true;

}

bool is_convergence(int index, vector<Chromosome> &chrome) {
	double sum[FUNC_NUM];
	for (int i = 0; i<FUNC_NUM; ++i)sum[i] = 0;
	int size = chrome.size();
	if (index == -1) {
		for (int i = 0; i<size - 1; ++i) {
			if (chrome[i] == chrome[i + 1]) {
			}
			else return false;
		}
	}
	else {
		for (int j = 0; j<size - 1; ++j) sum[index] += fabs(chrome[j + 1].f[index] - chrome[j].f[index]);
	}
	for (int i = 0; i<FUNC_NUM; ++i)
		if (sum[i]>0.1) return false;
	return true;
}

/*
初始化解5个解，并对每个解进行相应目标的优化，如：
先产生第1个解，
然后，对第1个解进行f[0]方向上的优化；
最后更新存档EP和边界解集
最终：EP.size()长度不确定，extreme.size() == 5
将5个边界值的相关目标设置为无穷大，如：第1个边界值的对应目标是f[0]

*/
void init_pop() {

	//将5个边界值的相关目标设置为无穷大，如：第1个边界值的对应目标是f[0] = INF
	for (int i = 0; i<FUNC_NUM; ++i)
		extreme[i].f[i] = INF;
	Chromosome chrome;				//初始解
	double r_lamda[FUNC_NUM];		//标志对哪个目标进行优化
	for (int i = 0; i < N; ++i) {			//产生5个解，并对每个解进行相应目标的优化
		clean_chromosome(chrome);
		init_route(chrome);			//随机生成新解，生成5个解
		compute_f(chrome);			//计算解的5个目标
		for (int j = 0; j < FUNC_NUM; ++j)
			r_lamda[j] = 0;			//对r_lamda数组进行初始化
		r_lamda[i] = 1;
		//第1个解，对f[0]进行搜索，产生新解；第2个解，对f[1]进行优化，产生新解
		local_search(chrome, r_lamda, -1);				//分别对五个目标进行优化
														//		update_EP(chrome,-1);
	}
}

void check();

vector<Chromosome> temp_EP;

void get_tempEP() {
	temp_EP.clear();
	temp_EP = EP;
	for (int i = 0; i<FUNC_NUM; ++i) {
		vector<Chromosome>::iterator iter;
		for (iter = temp_EP.begin(); iter != temp_EP.end(); ++iter) {
			if (is_better(*iter, extreme[i]) || is_equal(*iter, extreme[i])) break;
		}
		if (iter == temp_EP.end()) {
			iter = temp_EP.begin();
			while (iter != temp_EP.end()) {
				if (is_better(extreme[i], *iter)) iter = temp_EP.erase(iter);
				else ++iter;
			}
			temp_EP.push_back(extreme[i]);
		}
	}

}





/*
简单的冒泡排序，对某前沿面上的个体根据某个目标函数值进行升序排序
m1:每个前沿面上个体的数目；
fpara1[][2]: fpara1[][0],个体在种群中的索引
fpara1[][1]，个体在某个目标上的函数值。
*/
void sort1(int m1)
{
	float temp, temp1;
	int i1, j1, k1;
	for (k1 = 0; k1 < m1 - 1; k1++)
	{
		for (i1 = k1 + 1; i1 < m1; i1++)
		{
			if (fpara1[k1][1] > fpara1[i1][1])
			{
				temp = fpara1[k1][1];
				temp1 = fpara1[k1][0];
				fpara1[k1][1] = fpara1[i1][1];
				fpara1[k1][0] = fpara1[i1][0];
				fpara1[i1][1] = temp;
				fpara1[i1][0] = temp1;
			}
		}


	}
	return;
}

/*========================================================================*/
/*
rank	: 1,...,globalpop.maxrank
求拥挤距离,globalpop.cub_len[a]
逐个前沿求拥挤距离
1.获取每个前沿面上的个体数目:m1
2.循环遍历每个目标
前沿上个体根据相应目标,进行冒泡排序(升序)
对排序后的结果求值:
i = 0 || m1 - 1
length[i][0] = fpara1[i][0];
length[i][1] = 100*max;
其它
length[i][0] = fpara1[i][0];
length[i][1] = fabs(fpara1[i+1][1]- fpara1[i-1][1])
3. 求得每个解的,拥挤距离:globalpop.cub_len[a] += length[i][1];

*/
void gshare()
{

	float length[2 * maxpop][2], max;
	int i, j, m1, a;
	float min, Diff;  // Added 18.08.2003
					  //每个前沿面上的个体数目
	m1 = EP.size();

	for (j = 0; j < FUNC_NUM; j++)
	{
		for (i = 0; i < m1; i++)
		{
			fpara1[i][0] = 0;
			fpara1[i][1] = 0;
		}

		for (i = 0; i < m1; i++)
		{//种群中个体的索引
			fpara1[i][0] = i;
			//种群中第a个个体的第j个目标
			//fpara1[i][1] = globalpop.fitness[a][j];
			fpara1[i][1] = EP[i].f[j];
		}
		sort1(m1); /*Sort the arrays in ascending order of the fitness*/

		max = fpara1[m1 - 1][1];
		min = fpara1[0][1];  // Added 18.08.2003
		Diff = max - min;      // Added 18.08.2003 and 5 subsequent lines
		if (Diff < 0.0)
		{
			printf("Something wrong in keepaliven.h\n");
			exit(1);
		}
		else if (Diff == 0)
			Diff = 1;
		for (i = 0; i < m1; i++)
		{
			//对于边界解,设置无穷值
			if (i == 0 || i == (m1 - 1))
			{
				length[i][0] = fpara1[i][0];
				length[i][1] = 100 * max;
			}
			else	//对于中间解,等于相邻解的函数值差量的绝对规一化值
			{
				length[i][0] = fpara1[i][0];
				length[i][1] = fabs(fpara1[i + 1][1] - fpara1[i - 1][1]) / Diff; // crowding distances are normalized 18.08.2003
			}
		}
		for (i = 0; i < m1; i++)
		{
			a = (int)length[i][0];
			//求拥挤距离
			EP[a].cub_len += length[i][1];
		}
	}
	for (i = 0; i < m1; i++)
	{//种群中个体的索引
		fpara1[i][0] = i;
		//种群中第a个个体的第j个目标
		//fpara1[i][1] = globalpop.fitness[a][j];
		fpara1[i][1] = EP[i].cub_len;
	}
	sort1(m1); /*Sort the arrays in ascending order of the fitness*/
	vector<double> sp;
	for (i = 0; i < m1; i++)
	{
		EP[(int)fpara1[i][0]].ranking = m1 - i;

	}
	return;
}

void getMinMax()
{

	for (int i = 0; i < FUNC_NUM; i++)
	{
		double min = 999999, max = 0;
		for (int j = 0; j < EP.size(); j++)
		{
			if (EP[j].f[i] < min)
				min = EP[j].f[i];
			if (EP[j].f[i] > max)
				max = EP[j].f[i];
		}
		MMV[i][0] = min;
		if (min == max)
			MMV[i][1] = min + 1;
		else
			MMV[i][1] = max;

	}
}


void updateObjPossible()
{
	getMinMax();
	for (int i = 0; i < FUNC_NUM; i++)
	{
		if (MMV[i][1] - MMV[i][0] == 0)
			MMV[i][1] = 1 + MMV[i][0];
		//对某个目标的优化，会带来对其他目标的优化，增加了解的多样性，从而优化了其他目标
		for (int j = 0; j < EP.size(); j++)
		{
			FP[j][i] = 1 - ((EP[j].f[i] - MMV[i][0]) * 1.0) / (MMV[i][1] - MMV[i][0]);
		}
	}
}

/*
D:\3\running\movrptw_change\input\test50-0-0-0-0.d1.tw0DistanceMatrix.dat D:\3\running\movrptw_change\input\test50-0-0-0-0.d1.tw0TimeMatrix.dat D:\3\running\movrptw_change\input\test50-0-0-0-0.d1.tw0Specs.dat D:\3\running\movrptw_change\out\test_WheelSelection_i 3420 50

之前做了什么：
1.初始化所有客户点、客户点之间的行驶时间和行驶距离
2.拼凑了输出文件的路径，生成了结果输出文件的指针
*/

void process(double stopTime) {
	clock_t start, finish;
	start = clock();
	//初始化种群
	init_pop();

	Chromosome temp_chrome, chrome;
	double r_lamda[FUNC_NUM];
	int ind, ncount = 1;

	//记录每个局部搜索的次数
	//设定一个矩阵，记录第个LS对所有LS概率的影响
	for (int i = 0; i < FUNC_NUM; i++)
	{
		for (int j = 0; j < NEIGH_NUM; j++)				//领域操作
		{
			FIR[i][j] = 0;
			q[i][j] = 0;
		}
		LSArray[i][2] = 1.0 / 2.0;
		LSArray[i][1] = 1.0 / 3.0;
		LSArray[i][0] = 1.0 / 6.0;
	}

	for (int j = 0; j < NEIGH_NUM; j++)
	{
		LSTable[j] = false;
		lsn[j] = 0;
	}
	//更新个体的目标选择概率
	updateObjPossible();

	while (true) {

		clean_chromosome(temp_chrome);
		clean_chromosome(chrome);
		ind = rand() % EP.size();

		temp_chrome = EP[ind];
		chrome = temp_chrome;

		vector <double> sf;
		double sum = 0;
		for (int j = 0; j<FUNC_NUM; ++j)
		{
			r_lamda[j] = 0;
			//FP[ind][j]  += 0.000001;
			sum += FP[ind][j];
			sf.push_back(FP[ind][j]);
		}
		for (int j = 0; j < FUNC_NUM; j++)
			sf[j] /= sum;
		ind = RouletteWheelSelection(sf);
		r_lamda[ind] = 1;

		//更新LSArray表，以便局部搜索选择领域操作	


		local_search(temp_chrome, r_lamda, -1);

		int count = 0;
		double div;

		if (!is_equal(temp_chrome, chrome) && is_equal(temp_chrome, EP.back()))
			//更新成功了，一定会影响其中的解
		{	//更新最大值，最小值	
			//记录每个目标的改变量
			updateObjPossible();
		}


		finish = clock();
		if ((double)((finish - start) / CLOCKS_PER_SEC) > stopTime)
			break;
	}
}






void check(Chromosome &chrome) {
	bool used[MAX];
	vector<Route>::iterator r_iter;
	for (int j = 0; j <= cust_num; ++j) used[j] = false;
	for (r_iter = chrome.routes.begin(); r_iter != chrome.routes.end(); ++r_iter) {
		double t_time = 0;
		double w_time = 0;
		double capacity = 0;
		vector<int>::iterator c_iter;
		for (c_iter = r_iter->customers.begin(); c_iter != r_iter->customers.end() - 1; ++c_iter) {
			t_time += peer_time[*c_iter][*(c_iter + 1)];
			if (*(c_iter + 1) != 0) {
				if (t_time>customer[*(c_iter + 1)].e_time + max_delay_time) {
					printf("time false!!\n");
				}
			}
			else {
				if (t_time>customer[*(c_iter + 1)].e_time) {
					printf("end time false!!\n");
				}
			}
			w_time = t_time> customer[*(c_iter + 1)].b_time ? 0 : customer[*(c_iter + 1)].b_time - t_time;
			t_time += w_time + customer[*(c_iter + 1)].s_time;
			capacity += customer[*c_iter].demand;
			if (capacity>max_capacity) printf(" capacity false!\n");
			if (!used[*c_iter]) used[*c_iter] = true;
		}
	}
	int count = 0;
	for (int j = 1; j <= cust_num; ++j) if (used[j]) ++count;
	if (count != cust_num) {
		printf("less customer! \n");
	}
}

/*
对存档中的解进行排序，并依次输出每个解的5个目标值：
格式：
f[0] f[1] f[2] f[3] f[4]
......
*/
void output() {
	//给存档EP中的解，排序
	sort(EP.begin(), EP.end());

	//
	vector<Chromosome>::iterator iter;
	//for(iter=EP.begin(); iter!=EP.end();++iter){
	//	vector<Route>::iterator r_iter;
	//	for(r_iter=iter->routes.begin(); r_iter!=iter->routes.end();++r_iter){
	//		vector<int>::iterator c_iter;
	//		for(c_iter=r_iter->customers.begin();c_iter!=r_iter->customers.end();++c_iter){
	//			printf("%d ",*c_iter);
	//		}
	//		printf("\n");
	//	}
	//	printf("the route number is : %lf\n",iter->f[0]);
	//	printf("total distance is:%lf\n", iter->f[1]);
	//	printf("makespan is:%lf\n", iter->f[2]);
	//	printf("total waiting time is:%lf\n", iter->f[3]);
	//	printf("total delay time is:%lf\n\n",iter->f[4]);
	//}
	//printf("total solutions: %d\n", EP.size());
	cout << "输出：EP中解的" << EP.size() << "个目标值：" << endl;
	for (iter = EP.begin(); iter != EP.end(); ++iter) {
		for (int i = 0; i<FUNC_NUM; ++i) {
			printf("%lf ", iter->f[i]);
		}
		printf("\n");
	}

}

void output_best() {

	sort(total_best.begin(), total_best.end());
	vector<Chromosome>::iterator iter;
	for (iter = total_best.begin(); iter != total_best.end(); ++iter) {
		printf("%lf %lf %lf %lf %lf\n", iter->f[0], iter->f[1], iter->f[2], iter->f[3], iter->f[4]);
	}
	//	printf("total solutions: %d\n", total_best.size());

}

/*
对EP中的解进行检查，检查是不是正确解，是否是占优解
*/
void check() {

	bool used[MAX];
	//对每个解的路径进行检查，检查其客户点数是否等于cust_num 
	for (int i = 0; i < EP.size(); ++i) {													//遍历EP每个解
		vector<Route>::iterator r_iter;
		for (int j = 0; j <= cust_num; ++j)												//将客户点的标志数组初始化为false
			used[j] = false;
		//对解中路径进行时间约束和容量约束检查
		for (r_iter = EP[i].routes.begin(); r_iter != EP[i].routes.end(); ++r_iter) {		//遍历每个解中的路径
			double t_time = 0;
			double w_time = 0;
			double capacity = 0;
			vector<int>::iterator c_iter;
			for (c_iter = r_iter->customers.begin(); c_iter != r_iter->customers.end() - 1; ++c_iter) { //遍历每个路径中的客户点
				t_time += peer_time[*c_iter][*(c_iter + 1)];
				//路径r_iter上的时间约束检查
				if (*(c_iter + 1) != 0) {
					if (t_time > customer[*(c_iter + 1)].e_time + max_delay_time)
						printf("time false!!\n");
				}
				else {
					if (t_time > customer[*(c_iter + 1)].e_time)
						printf("end time false!!\n");
				}//时间约束检查通过
				w_time = t_time> customer[*(c_iter + 1)].b_time ? 0 : customer[*(c_iter + 1)].b_time - t_time;
				t_time += w_time + customer[*(c_iter + 1)].s_time;
				capacity += customer[*c_iter].demand;
				//容量约束检查
				if (capacity > max_capacity)
					printf(" capacity false!\n");
				if (!used[*c_iter])					//将遍历的客户点进行标记，每个解中只能用一次
					used[*c_iter] = true;
			}
		}
		int count = 0;
		for (int j = 1; j <= cust_num; ++j)
			if (used[j]) ++count;
		if (count != cust_num)						//如果遍历到的客户点数不等于总的客户点数
			printf("less customer! \n");
		//等于总的客户点数

		//printf("After check\n");
		//double dist=0;
		//for(r_iter=total_best[i].routes.begin(); r_iter!= total_best[i].routes.end(); ++r_iter){
		//	vector<int>::iterator c_iter;
		//	for(c_iter=r_iter->customers.begin(); c_iter!=r_iter->customers.end()-1; ++c_iter) dist+= peer_distance[*c_iter][*(c_iter+1)];
		//}
		//printf("distance is %lf\n", dist);
	}
	//对EP中的解检查其是否是占优解
	for (int i = 0; i<EP.size(); ++i)
		for (int j = 0; j<EP.size(); ++j)
			if (is_better(EP[i], EP[j])) {
				printf("%d is better than %d \n", i, j);
			}

}

//int main(){
//		char r_distance[100], r_time[100], r_spec[100];
//	//strcpy(r_distance,argv[1]);
//	//strcpy(r_time,argv[2]);
//	//strcpy(r_spec,argv[3]);
//	strcpy(r_distance,"distance.txt");
//	strcpy(r_time,"time.txt");
//	strcpy(r_spec,"spec.txt");
//	getData(r_distance,r_time,r_spec);
//	Route route;
//	route.customers.push_back(0);
//	route.customers.push_back(185);
//	route.customers.push_back(178);
//	route.customers.push_back(182);
//	route.customers.push_back(70);
////	route.customers.push_back(93);
//	route.customers.push_back(216);
//	route.customers.push_back(199);
//	route.customers.push_back(0);
//	compute_route_time(route);
//	 
//}

void initReadFile(char* c)
{


	vector<string> type;
	type.push_back("DistanceMatrix.dat");
	type.push_back("TimeMatrix.dat");
	type.push_back("Specs.dat");


	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{

			for (int k = 0; k < 5; k++)
			{

				for (int l = 0; l < 3; l++)
				{
					char path[100];
					char num[4];

					strcpy(path, c);
					strcat(path, "\\test");
					itoa(CP[i], num, 10);
					strcat(path, num);			//test50

					strcat(path, "-0-0-0-0.d");
					itoa(j, num, 10);
					strcat(path, num);			//test50-0-0-0-0.d0

					strcat(path, ".tw");
					itoa(k, num, 10);
					strcat(path, num);			//test50-0-0-0-0.d0.tw0

					strcat(path, type[l].c_str());
					files.push_back(path);
				}

			}
		}
	}
	return;
}

void readtime(char *c)
{
	fstream r_time(c, ios::in);

	for (int i = 0; i < 45; i++)
		r_time >> times[i];
}

/*
argc:代表参数的数量
控制台环境中程序名称后键入含参数的指令, 那么随后的参数将传递给argv[1]
输入数据：
argv:二维字符数组；
argv[1]:距离文件的路径字符数组
argv[2]:时间文件的路径字符数组
argv[3]:种群文件的路径字符数组
argv[4]：输出文件的根路径地址
argv[5]:stopTime字符串

D:\3\running\movrptw_change\input\test50-0-0-0-0.d1.tw0DistanceMatrix.dat D:\3\running\movrptw_change\input\test50-0-0-0-0.d1.tw0TimeMatrix.dat D:\3\running\movrptw_change\input\test50-0-0-0-0.d1.tw0Specs.dat D:\3\running\movrptw_change\out\5010_321\ 3420 50
D:\running\code\expercienceresult\mdls\input\test50-0-0-0-0.d1.tw0DistanceMatrix.dat D:\running\code\expercienceresult\mdls\input\test50-0-0-0-0.d1.tw0TimeMatrix.dat D:\running\code\expercienceresult\mdls\input\test50-0-0-0-0.d1.tw0Specs.dat D:\running\code\expercienceresult\mdls\out\1_321 2700 50
D:\running\code\expercienceresult\instance\Real-mdvrptw D:\running\code\expercienceresult\instance\times.txt D:\running\code\expercienceresult\instance\result\mdls_change2
*/

int main(int argc, char **argv) {				//一个int类型，一个二级指针。
	srand(time(NULL));							//得到机器的日历时间或者设置日历时间 
												//	 freopen("result.txt","w",stdout);
												//argc  >= 6
	clock_t start, finish;			//typedef long clock_t; long类型
	char r_distance[100], r_time[100], r_spec[100];		//定义char类型的三个一维数组		分别表示：	距离、时间、空间
														//strcpy(r_distance,argv[0]);
	initReadFile(argv[1]);			//D:\3\experimentalresult\instance\Real-mdvrptw,下面是各个.txt文件

	readtime(argv[2]);										//运行时间
	for (int n = 0; n < 45; n++)
	{
		strcpy(r_distance, files[n].c_str());								//argv[1]:距离文件的路径
		strcpy(r_time, files[n + 1].c_str());									//argv[2]:时间文件的路径
		strcpy(r_spec, files[n + 2].c_str());									//argv[3]:种群文件的路径
																				//atof()会扫描参数nptr字符串，跳过前面的空格字符，直到遇到数字或正负符号才开始做转换，
																				//而再遇到非数字或字符串结束时（'\0'）才结束转换，并将结果返回，str字符串可包含正负号、小数点或E(e)表示指数部分
		double stopTime = times[n];						   //把字符串转化成浮点数
														   //double stopTime=11400;
														   //功能:1.初始化所有客户点、客户点之间的行驶时间和行驶距离

		char out_file[100];
		strcpy(out_file, argv[3]);				//输出文件目录
		char num[4];
		itoa(n, num, 10);
		strcat(out_file, "\\");
		strcat(out_file, num);
		if (_access(out_file, 0) == -1)
		{
			_mkdir(out_file);
		}

		cout << "进行测试，输出r_distance、r_time、r_spec" << endl;
		cout << r_distance << endl;
		cout << r_time << endl;
		cout << r_spec << endl;

		cust_num = CP[n / 15];						//客户点数目 

		getData(r_distance, r_time, r_spec);
		int i = 0;
		start = clock();								//返回处理器调用某个进程或函数所花费的时间。
														//int i;
		for (i = 0; i < 30; ++i) {							//30表示:表示每个算例运行次数为30次

			char num[3];
			//功能：把一整数转换为字符串
			/*
			itoa(i,num,10);
			i ----需要转换成字符串的数字
			num---- 转换后保存字符串的变量
			10---- 转换数字的基数（即进制）。10就是说按10进制转换数字。还可以是2，8，16等等你喜欢的进制类型
			*/
			char out_file1[100];
			itoa(i, num, 10);
			strcpy(out_file1, out_file);
			strcat(out_file1, "\\");
			strcat(out_file1, num);			//outfile + "result" + i
			strcat(out_file1, ".txt");		//outfile + "result" + i + ".txt"
											/*
											freopen("in.txt","r",stdin);     //从in.txt 中读入数据
											freopen("out.txt","w",stdout);  // 将最后数据写入out.txt中 		*/
			freopen(out_file1, "w", stdout);		//以只写w方式打开文件out_file
													//在指定的时间内运行程序	

			process(stopTime);
			//对存档中的解进行排序，并依次输出每个解的5个目标值：	
			output();
			for (int j = 0; j < EP.size(); ++j)
				//功能：判断new_chromosome是否是占优解，删除其中相对new_chromosome不占优的，将占优的new_chromosome解加入占优解集total_best中
				update_best(EP[j]);
			//对EP中的解进行检查，检查是不是正确解，是否是占优解
			check();
			//EP解集清空
			EP.clear();
		}

		//将每个算例运行30次的总的最占优解集total_best输出到"final.txt"文件中
		finish = clock();
		char out_file1[100];
		strcpy(out_file1, out_file);
		strcat(out_file1, "\\");
		strcat(out_file1, "final.txt");
		freopen(out_file1, "w", stdout);
		output_best();
		//这个有必要吗
		check();
	}//n结束，第几个案例

	 //	printf("time: %lf\n",(double)((finish-start)/CLOCKS_PER_SEC)/i);
	 //cout << "结束了!" << endl;


	return 0;


}

