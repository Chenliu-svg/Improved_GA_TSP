/*
The GAlib-based genetic algorithm code for the Travelling Salesman Problem (TSP) .
作者：wying and Liuc
单位：华南理工大学软件学院
*/

#include <math.h>
#include <ctime>
#include "ga/GASStateGA.h"
#include "ga/GASimpleGA.h"
#include "ga/GA1DArrayGenome.h"
#include "ga/garandom.h"
#include "ga/std_stream.h"
#include<algorithm>
#include<set>
#include <cassert>
#include <string>
#define cout STD_COUT
#define cerr STD_CERR
#define endl STD_ENDL
#define ostream STD_OSTREAM
#define ifstream STD_IFSTREAM

// Set this up for your favorite TSP.  The sample one is a contrived problem
// with the towns laid out in a grid (so it is easy to figure out what the 
// shortest distance is, and there are many different paths with the same
// shortest path).  File format is that used by the TSPLIB problems.  You can 
// grab more problems from TSPLIB.
// 
#define MAX_TOWNS 100
#define TSP_FILE "./dataset/pr76.txt"
//#define TSP_FILE "berlin52.txt"
//#define TSP_FILE "lin105.txt"
//#define TSP_FILE "rat99.txt"

#define OUT_DIR "./outdir"

int x[MAX_TOWNS],y[MAX_TOWNS];//每个城市的x坐标和y坐标
int DISTANCE[MAX_TOWNS][MAX_TOWNS];//每两个城市之间的旅行成本，是对称的

// 算子申明
float TSPObjective(GAGenome&); //计算染色体的旅行总费用的目标函数
void  TSPInitializer(GAGenome&); //TSP问题的染色体初始化算子
int   Tranditional_TSPMutator(GAGenome&, float);//针对TSP问题的染色体变异算子
int   Favorable_TSPMutator(GAGenome&, float);//针对TSP问题的染色体变异算子
int   TSPCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);//针对TSP问题的染色体交叉算子

void writeTSPPath(ostream & os, GAGenome& g);//将指定染色体的旅行路线输出到指定文件


int try_mutate_times = 100;
int consistence_times = 60;
int main() {
  cout << "The GAlib program for the Travelling Salesman Problem (TSP) Berlin52.\n" << endl;


  //从Berlin52.txt文件读出各城市坐标
  double CityID;
  ifstream in(TSP_FILE); 
  if(!in) {
    cerr << "could not read data file " << TSP_FILE << "\n";
    exit(1);
  }
  int ntowns=0;
  do {
    in >> CityID;
    in >> x[ntowns];
    in >> y[ntowns];
    ntowns++;
  } while(!in.eof() && ntowns < MAX_TOWNS);
  in.close();
  if(ntowns >= MAX_TOWNS) {
    cerr << "data file contains more towns than allowed for in the fixed\n";
    cerr << "arrays.  Recompile the program with larger arrays or try a\n";
    cerr << "smaller problem.\n";
    exit(1);
  }

  //计算任意两个城市间的旅行成本
  double dx,dy;
  for(int i=0;i<ntowns;i++) {
    for(int j=i; j<ntowns;j++) {
      dx=x[i]-x[j]; dy=y[i]-y[j];
	  //注意取四舍五入之后的整数值
      DISTANCE[j][i]=DISTANCE[i][j]=floor(0.5 + sqrt(dx*dx+dy*dy) );
    }
  }

  //定义TSP问题的编码方案为一维的整数数组，其固定长度为城市个数
  GA1DArrayGenome<int> genome(ntowns);

  genome.evaluator(::TSPObjective);//为染色体指定计算目标值的函数
  genome.initializer(::TSPInitializer);//为染色体指定自定义的初始化算子
  genome.crossover(::TSPCrossover);//为染色体指定自定义的交叉算子
  genome.mutator(::Favorable_TSPMutator);//为染色体指定自定义的变异算子

  GASteadyStateGA ga(genome); 
  ga.nReplacement(16); 
  //选用稳态遗传算法进行TSP问题求解，指定其染色体编码方式
  //每一代要替换的个体数 = 16、总的运行代数500000，那么搜索的总个体数 = 16 * 500000 = 8000000
  ga.nGenerations(500000);
  ga.minimize();//为遗传算法指定优化目的是将目标函数值最小化
  ga.populationSize(200);//为遗传算法指定种群大小为200
  ga.pMutation(0.5);//为遗传算法指定变异概率
  ga.pCrossover(0.8);//为遗传算法指定交叉概率

  
  int epoch = 1;
  int best_epoch = 1;
  int best_res = 100000000;
  int worst_res = 0;
  int worst_epoch = 1;
  for (; epoch <= 10; epoch++) {
	  cout << "initializing..." << "\n"; cout.flush();
	  unsigned int seed = clock();
	  ga.initialize(seed);//使用从时钟得到的随机种子初始化遗传算法

	  cout <<"epoch:"<< epoch<< "  evolving..." << "\n"; cout.flush();
	  std::fstream fgacurve;
	  std::string curvepath= std::string(OUT_DIR)+ "/tspgacurve"+ std::to_string(epoch) +".txt";
	  fgacurve.open(curvepath, std::ios::out);
	  

	  int consistence = 0;
	  int cur_best = 100000000;
	  

	  //遗传算法开始迭代进化，直到达到指定的代数
	  while (!ga.done()) {
		  
		  ga.step();//进化一代
		  if (ga.generation() % (ga.nGenerations() / 50000) == 0)
		  {//进化过程中取100个采样点，记录进化过程中的最优目标值收敛信息到文件
			  int bestscore = ga.statistics().bestIndividual().score();
			  if (bestscore >= cur_best) {
				  consistence += 1;
			  }
			  else {
				  cur_best = bestscore;
				  consistence = 0;
			  }

			  cout << ga.generation() << "    " << bestscore << "\n"; cout.flush();
			  fgacurve << ga.generation() << "    " << bestscore << "\n";
		  }
		  // 如果连续300次没有进化，则退出 early stop
		  if (consistence >= consistence_times) {
			  cout << "连续300次进化没有提升，退出进化";
			  break;
		  }

	  }
	  int pre_score= ga.statistics().bestIndividual().score();
	  if (pre_score > worst_res) {
		  worst_res = pre_score;
		  worst_epoch = epoch;
	  }
	  if (pre_score < best_res) {
		  best_res = pre_score;
		  best_epoch = epoch;
	  }

	  fgacurve.close();

	  //遗传算法迭代终止后输出找到的最优旅行路线到文件
	  genome = ga.statistics().bestIndividual();
	  //cout << "\n" << "the shortest path found is "  << "\n";
	  //writeTSPPath(cout, genome);
	  std::fstream fbestpath;
	  std::string tsppathname = std::string(OUT_DIR) + "/tsppath" + std::to_string(epoch) + ".txt";
	  fbestpath.open(tsppathname, std::ios::out);
	  writeTSPPath(fbestpath, genome);
	  fbestpath.close();
	  cout << "the distance of the shortest path found: " << genome.score() << "\n";
  }

  cout << "best epoch: " << best_epoch << endl;
  cout << "worst epoch: " << worst_epoch << endl;
  return 0;
}


// Here are the genome operators that we want to use for this problem.
//计算染色体的旅行总费用的目标函数
float TSPObjective(GAGenome& g) {
  GA1DArrayGenome<int> & genome = (GA1DArrayGenome<int> &)g;
  int genomelength = genome.size();//genome.size()获取染色体的长度
  float dist = 0;
  int xx;
  int yy;

    for(int i=0; i<genomelength; i++) {
      xx = genome.gene(i);
	  // +1是因为到达最后一个城市后，返回第一个城市
      yy = genome.gene( (i+1)%genomelength );
	  // -1是因为数组下标从0开始
      dist += DISTANCE[xx-1][yy-1];
    }

  return dist;
}

//TSP问题的染色体初始化算子
void TSPInitializer(GAGenome& g) {
  GA1DArrayGenome<int> &genome=(GA1DArrayGenome<int> &)g;

  int genomelength = genome.size();
  int i,town;
  static bool visit[MAX_TOWNS];

  memset(visit, false, MAX_TOWNS*sizeof(bool));
  //GARandomInt(1,genomelength)生成1到genomelength之间的一个均匀随机整数
  town=GARandomInt(1,genomelength);
  // 随机取一个数，作为第一站
  visit[town-1]=true;
  genome.gene(0, town);//genome.gene(0, town)设置该染色体第0个基因位上的基因值为town
 
  for( i=1; i<genomelength; i++) {
    do {
		// 如果visit过了就重新选，相当于shuffle
      town=GARandomInt(1,genomelength);
    } while (visit[town-1]);
    visit[town-1]=true;
    genome.gene(i, town);
  }	
}
int Tranditional_TSPMutator(GAGenome& g, float pmut) {

	// 计算初始的适应度函数
	//float init_dist = TSPObjective(g);

	GA1DArrayGenome<int>& genome = (GA1DArrayGenome<int> &)g;
	int i;

	int genomelength = genome.size();
	float nmutator = pmut * genomelength;//要改变的边的数量

	int imutator = 0;
	while (imutator < nmutator) {
		/*float end_dist = 0;*/
		if (GARandomFloat() < 0.5) {//GARandomFloat()生成0到1之间的一个均匀随机浮点数
			//以0.5概率使用相互交换变异
			  //指在一个个体的染色体中随机选择两个位置，然后将这两个位置的基因进行互换

			
				int swapIndex1 = GARandomInt(0, genomelength - 1);
				int swapIndex2 = GARandomInt(0, genomelength - 1);
				int tmp;
				tmp = genome.gene(swapIndex2);
				genome.gene(swapIndex2, genome.gene(swapIndex1));
				genome.gene(swapIndex1, tmp);// swap only one time
			
			//计算交换以后的适应度函数，只有交换以后变好了，才允许变异


			imutator += 4;
		}
		else
		{
			//以0.5概率使用反转变异
			
				int inversion_start, inversion_end, tmp;
				inversion_start = GARandomInt(0, genomelength - 1);
				inversion_end = GARandomInt(0, genomelength - 1);
				// 保证start < end
				if (inversion_start > inversion_end)
				{
					tmp = inversion_start;
					inversion_start = inversion_end;
					inversion_end = tmp;
				}

				// 反转start 到 end的基因实现反转
				for (i = inversion_start; inversion_start < inversion_end; inversion_start++, inversion_end--)
				{
					// 
					tmp = genome.gene(inversion_start);
					genome.gene(inversion_start, genome.gene(inversion_end));
					genome.gene(inversion_end, tmp);
				}
			
			imutator += 2;
		}
	}

	return (1);
}




//针对TSP问题的染色体变异算子，pmut为变异概率
int Favorable_TSPMutator(GAGenome& g, float pmut) {
	
	// 计算初始的适应度函数
	float init_dist = TSPObjective(g);
	GA1DArrayGenome<int> & real = (GA1DArrayGenome<int> &)g;
	//GA1DArrayGenome<int> genome=(GA1DArrayGenome<int> &)g;
  int i;
  // 可以尝试try max time 次变异，找出最有利的变异：
  // 变异规则如下： 1. 好过父类：则直接选择该变异
  // 2. try max time次没有好于父类，则选择最优变异 
  // 从收敛代数和效果衡量
  int genomelength = real.size();
  float nmutator = pmut*genomelength;//要改变的边的数量

  int imutator=0;
  while( imutator<nmutator){
	  /*float end_dist = 0;*/
    if (GARandomFloat()<0.5) {//GARandomFloat()生成0到1之间的一个均匀随机浮点数
	  //以0.5概率使用相互交换变异
		//指在一个个体的染色体中随机选择两个位置，然后将这两个位置的基因进行互换
		GA1DArrayGenome<int> genome = real;
		GA1DArrayGenome<int> cur_best_genome = real;
		float cur_score = 100000000;
		float best_mutator_score = 100000000;
		int try_times = 0;
		do {
			genome = real;
			int swapIndex1 = GARandomInt(0, genomelength - 1);
			int swapIndex2 = GARandomInt(0, genomelength - 1);
			int tmp;
			tmp = genome.gene(swapIndex2);
			genome.gene(swapIndex2, genome.gene(swapIndex1));
			genome.gene(swapIndex1, tmp);// swap only one time
			cur_score= TSPObjective(genome);
			try_times = try_times + 1;
			
			if (try_times < try_mutate_times){
				if (best_mutator_score > cur_score) {
					best_mutator_score = cur_score;
					cur_best_genome = genome;
				}

			}
			else { genome = cur_best_genome; break; }
		

		} while (cur_score>=init_dist);
	  //计算交换以后的适应度函数，只有交换以后变好了，才允许变异
		real = genome;
	  
	  imutator+=4;
    }else
	  {
	  //以0.5概率使用反转变异
		
		GA1DArrayGenome<int> genome = real;
		GA1DArrayGenome<int> cur_best_genome = real;
		float cur_score = 100000000;
		float best_mutator_score = 100000000;
		int try_times = 0;
		do {
			genome = real;
			int inversion_start, inversion_end, tmp;
			inversion_start = GARandomInt(0, genomelength - 1);
			inversion_end = GARandomInt(0, genomelength - 1);
			// 保证start < end
			if (inversion_start > inversion_end)
			{
				tmp = inversion_start;
				inversion_start = inversion_end;
				inversion_end = tmp;
			}

			// 反转start 到 end的基因实现反转
			for (i = inversion_start; inversion_start < inversion_end; inversion_start++, inversion_end--)
			{
				// 
				tmp = genome.gene(inversion_start);
				genome.gene(inversion_start, genome.gene(inversion_end));
				genome.gene(inversion_end, tmp);
			}
			cur_score = TSPObjective(genome);
			try_times = try_times + 1;
			if (try_times < try_mutate_times) {
				if (best_mutator_score > cur_score) {
					best_mutator_score = cur_score;
					cur_best_genome = genome;
				}

			}
			else { genome = cur_best_genome; break; }

		} while (cur_score >= init_dist);
		real = genome;
	  imutator+=2;
    }
  }

  return (1);
}


//针对TSP问题的染色体交叉算子
int TSPCrossover(const GAGenome& g1, const GAGenome& g2, GAGenome* c1, GAGenome* c2) {
  GA1DArrayGenome<int> parent1=(GA1DArrayGenome<int> &)g1;
  GA1DArrayGenome<int> parent2=(GA1DArrayGenome<int> &)g2;

  int genomelength = parent1.size();

  int nc=0;// 统计交叉的数量

  GA1DArrayGenome<int> &child1=(GA1DArrayGenome<int> &)*c1;
  GA1DArrayGenome<int> &child2=(GA1DArrayGenome<int> &)*c2;

  if(c1)  {child1 = parent2; nc++;}
  if(c2)  {child2 = parent1; nc++;}

  
  //此处添加代码实现自己的交叉算子
  //交换两段基因
  int cross_start, cross_end, tmp, cross_len;
  
  cross_start = GARandomInt(0, genomelength - 1);
  cross_end = GARandomInt(0, genomelength - 1);
  //cross_end = cross_start+2;
  // 保证start < end
  if (cross_start > cross_end)
  {
	  tmp = cross_start;
	  cross_start = cross_end;
	  cross_end = tmp;
  }
  cross_len = cross_end-cross_start+1;

   //记录父代要交换的基因
  if (c1) {
	  // 创建一个set来保存交换的genes
	  std::set<int> cross_genes;
	  for (int k = cross_start; k <= cross_end; k++) {
		  cross_genes.insert(child1.gene(k));
	  }

	  //双指针遍历，i负责遍历p1,j负责遍历c1的空位置
	  for (int i = 0, j = 0; i < genomelength&&j< genomelength; i++) {
		  if (j == cross_start)j += cross_len;//跳过交叉区域位置
		  // 将parent1里不属于交换片段的基因拷贝到child1
		  if (!cross_genes.count(parent1.gene(i))){
			  
			  //child1.copy(parent1, j, i, 1);
			  child1.gene(j, parent1.gene(i));
			  assert(child1.gene(j) == parent1.gene(i));
			  
			  j++;
		  }
	  }
	  
  }
  assert(child1.size() == child2.size());
  std::set<int> set_child1;
  //判断没有重复城市出现
  for (int i = 0; i < genomelength; i++) {
	  set_child1.insert(child1.gene(i));
  }
  //cout <<"set_child1:" << set_child1.size() << endl;
  assert(set_child1.size() == genomelength);

 

  if (c2) {
	  // 创建一个set来保存交换的genes
	  std::set<int> cross_genes;
	  for (int k = cross_start; k <= cross_end; k++) {
		  cross_genes.insert(child2.gene(k));
	  }

	  //双指针遍历，i负责遍历p1,j负责遍历c1的空位置
	  for (int i = 0, j = 0; i < genomelength && j < genomelength; i++) {
		  if (j == cross_start)j += cross_len;//跳过交叉区域位置
		  // 将parent1里不属于交换片段的基因拷贝到child1
		  if (!cross_genes.count(parent2.gene(i))) {

			  //child1.copy(parent1, j, i, 1);
			  child2.gene(j, parent2.gene(i));
			  assert(child2.gene(j) == parent2.gene(i));

			  j++;
		  }
	  }

  }
  
  std::set<int> set_child2;
  //判断没有重复城市出现
  for (int i = 0; i < genomelength; i++) {
	  set_child2.insert(child2.gene(i));
  }
  //cout << "set_child2:" << set_child2.size() << endl;
  assert(set_child2.size() == genomelength);
  
  return nc;
}

//将指定染色体的旅行路线输出到指定文件
void writeTSPPath(ostream & os, GAGenome& g) 
{
	GA1DArrayGenome<int> & genome = (GA1DArrayGenome<int> &)g;
	int genomelength = genome.size();
	for(int i=0; i<genomelength; i++)
	{
		int xx = genome.gene(i);
		os << xx <<"    " <<x[xx-1] << "      "<<y[xx-1] << "\n";
	}
}
