/*
The GAlib-based genetic algorithm code for the Travelling Salesman Problem (TSP) .
���ߣ�wying and Liuc
��λ����������ѧ���ѧԺ
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

int x[MAX_TOWNS],y[MAX_TOWNS];//ÿ�����е�x�����y����
int DISTANCE[MAX_TOWNS][MAX_TOWNS];//ÿ��������֮������гɱ����ǶԳƵ�

// ��������
float TSPObjective(GAGenome&); //����Ⱦɫ��������ܷ��õ�Ŀ�꺯��
void  TSPInitializer(GAGenome&); //TSP�����Ⱦɫ���ʼ������
int   Tranditional_TSPMutator(GAGenome&, float);//���TSP�����Ⱦɫ���������
int   Favorable_TSPMutator(GAGenome&, float);//���TSP�����Ⱦɫ���������
int   TSPCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);//���TSP�����Ⱦɫ�彻������

void writeTSPPath(ostream & os, GAGenome& g);//��ָ��Ⱦɫ�������·�������ָ���ļ�


int try_mutate_times = 100;
int consistence_times = 60;
int main() {
  cout << "The GAlib program for the Travelling Salesman Problem (TSP) Berlin52.\n" << endl;


  //��Berlin52.txt�ļ���������������
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

  //���������������м�����гɱ�
  double dx,dy;
  for(int i=0;i<ntowns;i++) {
    for(int j=i; j<ntowns;j++) {
      dx=x[i]-x[j]; dy=y[i]-y[j];
	  //ע��ȡ��������֮�������ֵ
      DISTANCE[j][i]=DISTANCE[i][j]=floor(0.5 + sqrt(dx*dx+dy*dy) );
    }
  }

  //����TSP����ı��뷽��Ϊһά���������飬��̶�����Ϊ���и���
  GA1DArrayGenome<int> genome(ntowns);

  genome.evaluator(::TSPObjective);//ΪȾɫ��ָ������Ŀ��ֵ�ĺ���
  genome.initializer(::TSPInitializer);//ΪȾɫ��ָ���Զ���ĳ�ʼ������
  genome.crossover(::TSPCrossover);//ΪȾɫ��ָ���Զ���Ľ�������
  genome.mutator(::Favorable_TSPMutator);//ΪȾɫ��ָ���Զ���ı�������

  GASteadyStateGA ga(genome); 
  ga.nReplacement(16); 
  //ѡ����̬�Ŵ��㷨����TSP������⣬ָ����Ⱦɫ����뷽ʽ
  //ÿһ��Ҫ�滻�ĸ����� = 16���ܵ����д���500000����ô�������ܸ����� = 16 * 500000 = 8000000
  ga.nGenerations(500000);
  ga.minimize();//Ϊ�Ŵ��㷨ָ���Ż�Ŀ���ǽ�Ŀ�꺯��ֵ��С��
  ga.populationSize(200);//Ϊ�Ŵ��㷨ָ����Ⱥ��СΪ200
  ga.pMutation(0.5);//Ϊ�Ŵ��㷨ָ���������
  ga.pCrossover(0.8);//Ϊ�Ŵ��㷨ָ���������

  
  int epoch = 1;
  int best_epoch = 1;
  int best_res = 100000000;
  int worst_res = 0;
  int worst_epoch = 1;
  for (; epoch <= 10; epoch++) {
	  cout << "initializing..." << "\n"; cout.flush();
	  unsigned int seed = clock();
	  ga.initialize(seed);//ʹ�ô�ʱ�ӵõ���������ӳ�ʼ���Ŵ��㷨

	  cout <<"epoch:"<< epoch<< "  evolving..." << "\n"; cout.flush();
	  std::fstream fgacurve;
	  std::string curvepath= std::string(OUT_DIR)+ "/tspgacurve"+ std::to_string(epoch) +".txt";
	  fgacurve.open(curvepath, std::ios::out);
	  

	  int consistence = 0;
	  int cur_best = 100000000;
	  

	  //�Ŵ��㷨��ʼ����������ֱ���ﵽָ���Ĵ���
	  while (!ga.done()) {
		  
		  ga.step();//����һ��
		  if (ga.generation() % (ga.nGenerations() / 50000) == 0)
		  {//����������ȡ100�������㣬��¼���������е�����Ŀ��ֵ������Ϣ���ļ�
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
		  // �������300��û�н��������˳� early stop
		  if (consistence >= consistence_times) {
			  cout << "����300�ν���û���������˳�����";
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

	  //�Ŵ��㷨������ֹ������ҵ�����������·�ߵ��ļ�
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
//����Ⱦɫ��������ܷ��õ�Ŀ�꺯��
float TSPObjective(GAGenome& g) {
  GA1DArrayGenome<int> & genome = (GA1DArrayGenome<int> &)g;
  int genomelength = genome.size();//genome.size()��ȡȾɫ��ĳ���
  float dist = 0;
  int xx;
  int yy;

    for(int i=0; i<genomelength; i++) {
      xx = genome.gene(i);
	  // +1����Ϊ�������һ�����к󣬷��ص�һ������
      yy = genome.gene( (i+1)%genomelength );
	  // -1����Ϊ�����±��0��ʼ
      dist += DISTANCE[xx-1][yy-1];
    }

  return dist;
}

//TSP�����Ⱦɫ���ʼ������
void TSPInitializer(GAGenome& g) {
  GA1DArrayGenome<int> &genome=(GA1DArrayGenome<int> &)g;

  int genomelength = genome.size();
  int i,town;
  static bool visit[MAX_TOWNS];

  memset(visit, false, MAX_TOWNS*sizeof(bool));
  //GARandomInt(1,genomelength)����1��genomelength֮���һ�������������
  town=GARandomInt(1,genomelength);
  // ���ȡһ��������Ϊ��һվ
  visit[town-1]=true;
  genome.gene(0, town);//genome.gene(0, town)���ø�Ⱦɫ���0������λ�ϵĻ���ֵΪtown
 
  for( i=1; i<genomelength; i++) {
    do {
		// ���visit���˾�����ѡ���൱��shuffle
      town=GARandomInt(1,genomelength);
    } while (visit[town-1]);
    visit[town-1]=true;
    genome.gene(i, town);
  }	
}
int Tranditional_TSPMutator(GAGenome& g, float pmut) {

	// �����ʼ����Ӧ�Ⱥ���
	//float init_dist = TSPObjective(g);

	GA1DArrayGenome<int>& genome = (GA1DArrayGenome<int> &)g;
	int i;

	int genomelength = genome.size();
	float nmutator = pmut * genomelength;//Ҫ�ı�ıߵ�����

	int imutator = 0;
	while (imutator < nmutator) {
		/*float end_dist = 0;*/
		if (GARandomFloat() < 0.5) {//GARandomFloat()����0��1֮���һ���������������
			//��0.5����ʹ���໥��������
			  //ָ��һ�������Ⱦɫ�������ѡ������λ�ã�Ȼ��������λ�õĻ�����л���

			
				int swapIndex1 = GARandomInt(0, genomelength - 1);
				int swapIndex2 = GARandomInt(0, genomelength - 1);
				int tmp;
				tmp = genome.gene(swapIndex2);
				genome.gene(swapIndex2, genome.gene(swapIndex1));
				genome.gene(swapIndex1, tmp);// swap only one time
			
			//���㽻���Ժ����Ӧ�Ⱥ�����ֻ�н����Ժ����ˣ����������


			imutator += 4;
		}
		else
		{
			//��0.5����ʹ�÷�ת����
			
				int inversion_start, inversion_end, tmp;
				inversion_start = GARandomInt(0, genomelength - 1);
				inversion_end = GARandomInt(0, genomelength - 1);
				// ��֤start < end
				if (inversion_start > inversion_end)
				{
					tmp = inversion_start;
					inversion_start = inversion_end;
					inversion_end = tmp;
				}

				// ��תstart �� end�Ļ���ʵ�ַ�ת
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




//���TSP�����Ⱦɫ��������ӣ�pmutΪ�������
int Favorable_TSPMutator(GAGenome& g, float pmut) {
	
	// �����ʼ����Ӧ�Ⱥ���
	float init_dist = TSPObjective(g);
	GA1DArrayGenome<int> & real = (GA1DArrayGenome<int> &)g;
	//GA1DArrayGenome<int> genome=(GA1DArrayGenome<int> &)g;
  int i;
  // ���Գ���try max time �α��죬�ҳ��������ı��죺
  // ����������£� 1. �ù����ࣺ��ֱ��ѡ��ñ���
  // 2. try max time��û�к��ڸ��࣬��ѡ�����ű��� 
  // ������������Ч������
  int genomelength = real.size();
  float nmutator = pmut*genomelength;//Ҫ�ı�ıߵ�����

  int imutator=0;
  while( imutator<nmutator){
	  /*float end_dist = 0;*/
    if (GARandomFloat()<0.5) {//GARandomFloat()����0��1֮���һ���������������
	  //��0.5����ʹ���໥��������
		//ָ��һ�������Ⱦɫ�������ѡ������λ�ã�Ȼ��������λ�õĻ�����л���
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
	  //���㽻���Ժ����Ӧ�Ⱥ�����ֻ�н����Ժ����ˣ����������
		real = genome;
	  
	  imutator+=4;
    }else
	  {
	  //��0.5����ʹ�÷�ת����
		
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
			// ��֤start < end
			if (inversion_start > inversion_end)
			{
				tmp = inversion_start;
				inversion_start = inversion_end;
				inversion_end = tmp;
			}

			// ��תstart �� end�Ļ���ʵ�ַ�ת
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


//���TSP�����Ⱦɫ�彻������
int TSPCrossover(const GAGenome& g1, const GAGenome& g2, GAGenome* c1, GAGenome* c2) {
  GA1DArrayGenome<int> parent1=(GA1DArrayGenome<int> &)g1;
  GA1DArrayGenome<int> parent2=(GA1DArrayGenome<int> &)g2;

  int genomelength = parent1.size();

  int nc=0;// ͳ�ƽ��������

  GA1DArrayGenome<int> &child1=(GA1DArrayGenome<int> &)*c1;
  GA1DArrayGenome<int> &child2=(GA1DArrayGenome<int> &)*c2;

  if(c1)  {child1 = parent2; nc++;}
  if(c2)  {child2 = parent1; nc++;}

  
  //�˴���Ӵ���ʵ���Լ��Ľ�������
  //�������λ���
  int cross_start, cross_end, tmp, cross_len;
  
  cross_start = GARandomInt(0, genomelength - 1);
  cross_end = GARandomInt(0, genomelength - 1);
  //cross_end = cross_start+2;
  // ��֤start < end
  if (cross_start > cross_end)
  {
	  tmp = cross_start;
	  cross_start = cross_end;
	  cross_end = tmp;
  }
  cross_len = cross_end-cross_start+1;

   //��¼����Ҫ�����Ļ���
  if (c1) {
	  // ����һ��set�����潻����genes
	  std::set<int> cross_genes;
	  for (int k = cross_start; k <= cross_end; k++) {
		  cross_genes.insert(child1.gene(k));
	  }

	  //˫ָ�������i�������p1,j�������c1�Ŀ�λ��
	  for (int i = 0, j = 0; i < genomelength&&j< genomelength; i++) {
		  if (j == cross_start)j += cross_len;//������������λ��
		  // ��parent1�ﲻ���ڽ���Ƭ�εĻ��򿽱���child1
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
  //�ж�û���ظ����г���
  for (int i = 0; i < genomelength; i++) {
	  set_child1.insert(child1.gene(i));
  }
  //cout <<"set_child1:" << set_child1.size() << endl;
  assert(set_child1.size() == genomelength);

 

  if (c2) {
	  // ����һ��set�����潻����genes
	  std::set<int> cross_genes;
	  for (int k = cross_start; k <= cross_end; k++) {
		  cross_genes.insert(child2.gene(k));
	  }

	  //˫ָ�������i�������p1,j�������c1�Ŀ�λ��
	  for (int i = 0, j = 0; i < genomelength && j < genomelength; i++) {
		  if (j == cross_start)j += cross_len;//������������λ��
		  // ��parent1�ﲻ���ڽ���Ƭ�εĻ��򿽱���child1
		  if (!cross_genes.count(parent2.gene(i))) {

			  //child1.copy(parent1, j, i, 1);
			  child2.gene(j, parent2.gene(i));
			  assert(child2.gene(j) == parent2.gene(i));

			  j++;
		  }
	  }

  }
  
  std::set<int> set_child2;
  //�ж�û���ظ����г���
  for (int i = 0; i < genomelength; i++) {
	  set_child2.insert(child2.gene(i));
  }
  //cout << "set_child2:" << set_child2.size() << endl;
  assert(set_child2.size() == genomelength);
  
  return nc;
}

//��ָ��Ⱦɫ�������·�������ָ���ļ�
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
