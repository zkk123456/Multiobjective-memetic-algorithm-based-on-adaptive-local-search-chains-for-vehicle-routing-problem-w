#include<vector>
using namespace std;

#ifndef VRPTW_H
#define VRPTW_H
#define maxpop   500  /*Max population */

const int MAX = 300;
const int FUNC_NUM = 5;
const int NEIGH_NUM = 3;
const int INF = 9999999;


//客户点
static struct Customer {
	int id;					//编号
	double x;				//x轴度坐标
	double y;				//y轴坐标
	double s_time;			//服务时间	
	double b_time;			//最早服务时间	
	double e_time;			//最晚服务时间
	double demand;			//需求量
}customer[MAX];				//客户点数组

double a_time[MAX];                         //arrival time arrays
double l_time[MAX];                         //leave time arrays
double w_time[MAX];                         //wait time	arrays
double d_time[MAX];                          //delay time	arrays
double total_wait[MAX];                     //the total wait time before i
double total_delay[MAX];                     //the total delay time before i

double max_wait[MAX];						//the max wait time最大合法延迟时间

double b_demand[MAX];                      // the total demand before customer i including i in the route
double a_demand[MAX];                      // the total demand after customer i including i	in the route

										   //求i点之前之后的路径长度有什么用
double b_distance[MAX];                    //the total distance before customer i including i
double a_distance[MAX];                    //the total distance after customer i including i


struct Route {							//结构体：路径
	vector<int> customers;				//路径中含有的客户点标号
	double travel_dist;					//这条路径的行驶距离
	double travel_time;					//这条路径的行驶时间
	double wait_time;					//这条路径的等待时间
	double delay_time;					//这条路径的延迟时间
	double capacity;					//容量
};

int cust_num;							//客户点数目
double max_capacity = 200;				//最大容量
double max_delay_time = 1800;				//最大延迟时间1800


double peer_time[MAX][MAX];				//两个客户点之间的行驶时间
double peer_distance[MAX][MAX];			//两个客户点之间的行驶距离

int archive_size = 500;					//存档大小.
										/*double epsilon[FUNC_NUM];  */            //epsilon 保存极值
double epsilon = 0.05;					//基于<-占优关系的存档参数：


const int N = 5;							//目标数目



struct Chromosome {						//解
	vector<Route> routes;				//路径集合
	double f[5];						//5个目标值,f[0]:车辆数目,f[1]:总行驶距离;f[2]:最长路径行驶时间;f[3]:总等待时间;f[4]:总延迟时间
	int box_f[5];						//非占优解的关联矩阵代表解在<-解空间上的位置
										//int s_f[FUNC_NUM];					//成功、失败的次数
	int local[FUNC_NUM];					//经过的局部搜索
	int count;
	double cub_len;
	int ranking;

	bool operator < (const Chromosome &chrome) {    //是要根据前三个目标来确定大小吗？
		return f[0]<chrome.f[0] || (f[0] == chrome.f[0] && f[1]<chrome.f[1]) || (f[0] == chrome.f[0] && f[1] == chrome.f[1] && f[2]<chrome.f[2]);
	}

	bool operator == (const Chromosome &chrome) {
		bool used[MAX][MAX];
		for (int i = 0; i <= cust_num; ++i)
			for (int j = 0; j <= cust_num; ++j)
				used[i][j] = false;

		for (vector<Route>::iterator r_iter = routes.begin(); r_iter != routes.end(); ++r_iter)
			for (vector<int>::iterator c_iter = r_iter->customers.begin(); c_iter != r_iter->customers.end() - 1; ++c_iter)
				used[*c_iter][*(c_iter + 1)] = !used[*c_iter][*(c_iter + 1)];					//将解中路径中相连接的客户点标志为true

		for (vector<Route>::const_iterator r_iter = chrome.routes.begin(); r_iter != chrome.routes.end(); ++r_iter)
			for (vector<int>::const_iterator c_iter = r_iter->customers.begin(); c_iter != r_iter->customers.end() - 1; ++c_iter)
				used[*c_iter][*(c_iter + 1)] = !used[*c_iter][*(c_iter + 1)];

		for (int i = 0; i <= cust_num; ++i)
			for (int j = 0; j <= cust_num; ++j)
				if (used[i][j]) return false;
		return true;
	}
	double similarity;									//类似; 相像性; 相仿性; 类似性，相似物;
	double crowding;									//拥挤现象; 晕线加密;
};


struct SUCCESS {							//结构体：路径
	vector<Chromosome> customers;				//路径中含有的客户点标号
												//double *sf[1000][FUNC_NUM];					//成功的次数比例
	int sf[500][FUNC_NUM];
};

vector<Chromosome> total_best;							//最好的解集合
Chromosome extreme[FUNC_NUM];							//边界解集合：代表5个目标上最优的解；如：extreme[0]:f[0]上最优，extreme[1]:f[1]上最优，extreme[2]:f[2]上最优，extreme[3]:f[3]上最优，extreme[4]:f[4]上最优
SUCCESS SIA;										//成功信息存档

vector<Chromosome>chromosome;						   //解集
vector<Chromosome> EP;								  //
													  //vector<Chromosome> BETWEEN;								//中间解
vector<double> IS;
double FP[maxpop][FUNC_NUM];
double MMV[FUNC_NUM][2];									//SFC[FUNC_NUM][0]:目标的最小值 ;SFC[FUNC_NUM][1]:目标的最大值
double LSArray[FUNC_NUM][NEIGH_NUM];
double FIR[FUNC_NUM][NEIGH_NUM], q[FUNC_NUM][NEIGH_NUM];
int lsn[NEIGH_NUM];

bool LSTable[NEIGH_NUM];
vector<Chromosome> parent;								//上一代解
vector<Chromosome> children;								//子代解
vector< vector<int> > F;								//放int值集合的集合
vector< vector<int>> S;									//放int值集合的集合

vector<string> files;

double times[45];
int ITER = 3000;											//子项目数3000
int ind2 = 0;											//记录选择的领域
int CP[3] = { 50,150,250 };
bool EP_flag, is_EP = false;											//EP标志

float fpara1[2 * maxpop][2];

#endif