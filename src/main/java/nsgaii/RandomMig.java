package nsgaii;

import networkInit.DeviceInit;
//import networkInit.ServiceInit;
import networkInit.Round;
import networkInit.ServiceInit;
import variable.Variable;

import java.util.*;
//import java.util.Vector;


//import variable.Variable;
//import java.util.Scanner;

public class RandomMig {
	private int maxGen; //最大进化代数
		private int gen;// 当前进化代数

	private int lowBand;
	private int highBand;
	double pc;// 交叉概率
	double pm;// 变异概率
	
	ArrayList<Chromsome> popList;
	ArrayList<Chromsome> prePopList;
	ArrayList<ServiceInit> service;
	ArrayList<DeviceInit> device;
	Vector<Integer> serviceReq;
	Round iteRound;
	int numService;
	int popSize;
	
	Random random = new Random();
	ArrayList<Chromsome> randomPopList;

	private int[][] migArray ;// 放置迁移矩阵
	private int[] migration;// 放置染色体数组
	HashMap<Chromsome,ArrayList<Chromsome>> dominatingMap;//支配集
	
	public RandomMig(int popSize,int maxGen,double pc, double pm,ArrayList<Chromsome> popList,int numService,ArrayList<ServiceInit> service,ArrayList<DeviceInit> device,Vector<Integer> serviceReq,Round iteRound,int lowBand,int highBand) {
		super();
		this.popSize = popSize;
		this.maxGen = maxGen;
		this.pc = pc;
		this.pm = pm;
		this.popList = popList;
		this.numService = numService;
		this.service = service;
		this.device = device;
		this.serviceReq = serviceReq;
		this.iteRound = iteRound;
		this.lowBand = lowBand;
		this.highBand = highBand;
	

		this.prePopList = new ArrayList<Chromsome>();
		for (int i = 0; i < this.popSize; i++) {
			Chromsome chromsome = new Chromsome(7);
			prePopList.add(chromsome);
		}	
	}
	/*
	 * 随机迁移设备算法
	 */
	// 判断容器的数量是否超过了5个
	public boolean JudegIsSuit(int[][] crr) {
		int sum = 0;
		for (int i = 0; i < crr.length; i++) {
			for (int j = 0; j < crr[i].length; j++) {
				sum += crr[i][j];
			}
			if (sum > 5) {
				return true;
			}
			sum = 0;
		}
		return false;
	}
	// 随机生成迁移矩阵
	public int[][] RandomArray() {
		int Mig[][] = new int[numService][numService];// 初始化二维数组，此时的二维数组设为numService
		Random random = new Random();
		for (int i = 0; i < numService; ++i) {
			int a = random.nextInt(numService);// 在numService内随机生成一个数
			for (int j = 0; j < numService; ++j) {
				// 当j值相等时，赋值给本列的某一行
				if (j == a)
					Mig[j][i] = 1;
				else
					Mig[j][i] = 0;
			}
		}
		// 输出生成的随机矩阵
//		for (int i = 0; i < numService; ++i) {
//			for (int j = 0; j < numService; ++j) {
//				System.out.print(Mig[i][j] + " ");
//			}
//			System.out.println();
//		}
//		System.out.println("完成矩阵的产生");
		if (JudegIsSuit(Mig) == true) // 如果不满足限制条件，重新生成随机矩阵
			Mig = RandomArray();
		return Mig;
	}
	// 矩阵乘法得到新的迁移数组
	public int[] NewArray(ArrayList<Integer> x, int[][] arr) {
		int[] Final = new int[numService];
		for (int i = 0; i < numService; i++) {
			for (int j = 0; j < numService; j++) {
				Final[i] += x.get(j) * arr[j][i];
			}
		}
		return Final;
	}
	// 计算两点之间的欧氏距离
	public float calDij(DeviceInit device1, DeviceInit device2) {
		return (float) Math.sqrt((device1.getX() - device2.getX()) * (device1.getX() - device2.getX())
				+ (device1.getY() - device2.getY()) * (device1.getY() - device2.getY()));
	}

	// 计算两点之间的信噪比
	public float calNoise(DeviceInit device1, DeviceInit device2) {
		if (device1.getFlag() == 1 && device2.getFlag() == 1) {
			return (float) (Math.pow(10, ((-174 + 10 * (Math.log(lowBand) / Math.log(10))) / 10)) / 1000);
		} else
			return (float) (Math.pow(10, ((-174 + 10 * (Math.log(highBand) / Math.log(10))) / 10)) / 1000);
	}

	// 计算两个物理节点之间的rij
	public float calRij(DeviceInit device1, DeviceInit device2) {
		float rij = 0.0f;
		if (device1.getFlag() == 1 && device2.getFlag() == 1) {
			rij = (float) (lowBand
					* Math.log(1 + 0.2 * Math.pow(calDij(device1, device2), -4) / calNoise(device1, device2))
					/ Math.log(2));
		} else {
			rij = (float) (highBand
					* Math.log(1 + 0.2 * Math.pow(calDij(device1, device2), -4) / calNoise(device1, device2))
					/ Math.log(2));
		}
//		System.out.println("测试0的值：");
//		System.out.println(rij);
		return rij;
	}

	public void copyPopList() {
		//init.initPopulation();// 生成初始种群
		randomPopList = new ArrayList<Chromsome>();
		//复制PopList；
		for (int i = 0; i < popList.size(); i++) {
			Chromsome chromsome = new Chromsome(7);
			randomPopList.add(chromsome);
			//randomPopList.set(i, init.popList.get(i));
			randomPopList.set(i, (Chromsome)popList.get(i).clone());
		}
	}		
	// 初始化种群染色体，得到最终的经过迁移的染色体
	// 算法运行逻辑
	public void FinallyChormsome() {
		migArray = new int[numService][numService];
		migration = new int[numService];
		copyPopList();
		float c_cst = 0.0f;// 计算内存限制
		float e_trs = 0.0f;// 计算能耗限制
		// 输出没有经过限制条件的随机迁移情况
		for (int i = 0; i < popList.size(); i++) {
			migArray = RandomArray();
			migration = NewArray(popList.get(i).deviceId, migArray).clone();
			for (int j = 0; j < numService; j++) {
				randomPopList.get(i).deviceId.set(j, migration[j]);
			}		
		}
		// 输出染色体的种群
		//System.out.println("初始生成的随机迁移情况：");
		/*for (int i = 0; i < popList.size(); i++) {
			for (int j = 0; j < numService; j++) {
				System.out.print(randomPopList.get(i).deviceId.get(j) + " ");
			}
			System.out.println();
		}*/
		//System.out.println("经过限制条件以后的染色体：");
		for (int i = 0; i < popList.size(); i++) {
			for (int j = 0; j < numService; j++) {
				for (int m = 0; m < numService; m++) {
					if (randomPopList.get(i).deviceId.get(j).equals(randomPopList.get(i).deviceId.get(m))) 
						c_cst += service.get(serviceReq.get(m)).getWkd();
					if (popList.get(i).deviceId.get(j).equals(randomPopList.get(i).deviceId.get(m)) && j!=m){
							e_trs += 0.2 * service.get(serviceReq.get(m)).getWkd() * 1048576
									/ calRij(device.get(popList.get(i).deviceId.get(m)), 
											device.get(popList.get(i).deviceId.get(j)));
							//System.out.println("第一台设备是："+init.getDevice().get(init.getPopList().get(i).deviceId.get(m)).getId());
							//System.out.println("第二台设备是："+init.getDevice().get(init.getPopList().get(i).deviceId.get(j)).getId());
//							Systedevice = {ArrayList@511}  size = 20m.out.println("计算一下服务的能耗：" );
//							System.out.println(e_trs);
					}
				}
				//当某个设备的传输能耗以及服务大小之和大于最大剩余能量和最大剩余容量时
				while (c_cst > device.get(j).getStg()  || e_trs > device.get(j).getRsd()) {
					migArray = RandomArray();
					migration = NewArray(popList.get(i).deviceId, migArray).clone();
					//将新产生的矩阵赋值到染色体中
					for(int a=0;a<numService;a++) {
						randomPopList.get(i).deviceId.set(a, migration[a]);
					}
					//重新计算是否有超出限制的情况发生
					for (int p = 0; p < numService; p++) {			
						c_cst = 0.0f;
						e_trs = 0.0f;
						for (int q = 0; q < numService; q++) {
							if (popList.get(i).deviceId.get(p).equals(randomPopList.get(i).deviceId.get(q))) 
								c_cst += service.get(serviceReq.get(q)).getWkd();
							if (popList.get(i).deviceId.get(p).equals(randomPopList.get(i).deviceId.get(q)) && p!=q){
									e_trs += 0.2 * service.get(serviceReq.get(q)).getWkd() * 1048576
											/ calRij(device.get(popList.get(i).deviceId.get(q)), 
													device.get(popList.get(i).deviceId.get(p)));
							}
						}
						
					}
				}
			c_cst = 0.0f;
			e_trs = 0.0f;
			}	
		}

		for(int i = 0; i < randomPopList.size() ; i ++) {
			for(int j =0 ; j < randomPopList.get(i).deviceId.size(); j++) {
				popList.get(i).mgDeviceId.set(j,randomPopList.get(i).deviceId.get(j));
			}
		}
/*		System.out.println("-------------poplist * numservice-------------");
		System.out.println(popList.size());
		System.out.println(numService);*/
/*		for (int i = 0; i < popList.size(); i++) {
			for (int j = 0; j < numService; j++) {
				System.out.print(randomPopList.get(i).deviceId.get(j) + " ");
			}
			System.out.println();
		}*/
		//return randomPopList;
	}
	
	
	//计算迁移时间和迁移的能耗
	public void calinsideFitness() {
		FinallyChormsome();
		double t = 0.0;
		double e = 0.0;
		double cst = 0.0;
		double cbf = 0.0;
		for(int i = 0;i < popList.size();i++) {
				//double[] fitness = new double[3];
				for(int m = 0;m < numService;m++) {
					for(int n = 0;n < numService;n++) {
						if(popList.get(i).deviceId.get(m).equals(randomPopList.get(i).deviceId.get(n))) {
							cst += service.get(serviceReq.get(n)).getWkd();	
						}
					}
					cbf += (device.get(popList.get(i).deviceId.get(m)).getRsd() - cst) / device.get(popList.get(i).deviceId.get(m)).getRsd();
					cst = 0.0;
					if(popList.get(i).deviceId.get(m) != randomPopList.get(i).deviceId.get(m)) {
						t += service.get(serviceReq.get(m)).getWkd() * 1048576 / calRij(device.get(popList.get(i).deviceId.get(m)), device.get(randomPopList.get(i).deviceId.get(m)))
								+ service.get(serviceReq.get(m)).getCr() / device.get(randomPopList.get(i).deviceId.get(m)).getF();
						e += 0.2 * service.get(serviceReq.get(m)).getWkd() * 1048576/ calRij(device.get(popList.get(i).deviceId.get(m)), device.get(randomPopList.get(i).deviceId.get(m)))
								+ 0.002 * service.get(serviceReq.get(m)).getCr() / device.get(randomPopList.get(i).deviceId.get(m)).getF();
					}
					else {
						if(m < (numService - 1)) {
							t += service.get(serviceReq.get(m)).getCr() / device.get(randomPopList.get(i).deviceId.get(m)).getF()
									+ service.get(serviceReq.get(m)).getDt() / calRij(device.get(popList.get(i).deviceId.get(m)), device.get(popList.get(i).deviceId.get(m+1)));
							e += 0.2 * service.get(serviceReq.get(m)).getDt() / calRij(device.get(popList.get(i).deviceId.get(m)), device.get(popList.get(i).deviceId.get(m+1)))
									+ 0.5 * service.get(serviceReq.get(m)).getCr() / device.get(randomPopList.get(i).deviceId.get(m)).getF();
						}
						else if(m == (numService - 1)){
							t += service.get(serviceReq.get(m)).getCr() / device.get(randomPopList.get(i).deviceId.get(m)).getF();
							e += 0.5 * service.get(serviceReq.get(m)).getCr() / device.get(randomPopList.get(i).deviceId.get(m)).getF();
						}
						
					}
				}
				cbf = cbf / numService;
				popList.get(i).insidefitness[0] = t;
				popList.get(i).insidefitness[1] = e;
				popList.get(i).insidefitness[2] = cbf;
				//	System.out.print("t=" + fitness[0] + " " + "e=" + fitness[1] + " " + "cbf=" + fitness[2]);
				//System.out.println();
				t = 0.0;
				e = 0.0;
				cbf = 0.0;
			}
	}
	/*
	 * 以下是跑外层Nsga算法的流程
	 * 
	 * 以下是实验的结果
	 */
	//计算两个设备之间的空间距离
	public static float evaluateLoc(DeviceInit device, Round iteRound) {
		float locRet = 0;
		float d = (float) Math.sqrt((iteRound.getX() - device.getX()) * (iteRound.getX() - device.getX())
				+ (iteRound.getY() - device.getY()) * (iteRound.getY() - device.getY()));
		if (d <= Math.abs(iteRound.getRadius() - Variable.r)) {
			float r = iteRound.getRadius() < Variable.r ? iteRound.getRadius() : Variable.r;
			return (float) ((Math.PI * r * r) / iteRound.getArea());
		}
		float ang1 = (float) Math.acos((iteRound.getRadius() * iteRound.getRadius() + d * d - Variable.r * Variable.r)
				/ 2. / iteRound.getRadius() / d);
		float ang2 = (float) Math.acos(
				(Variable.r * Variable.r + d * d - iteRound.getRadius() * iteRound.getRadius()) / 2. / Variable.r / d);
		float ret = (float) (ang1 * iteRound.getRadius() * iteRound.getRadius() + ang2 * Variable.r * Variable.r
				- d * iteRound.getRadius() * Math.sin(ang1));
		locRet += (float) (ret / (Math.PI * Variable.r * Variable.r));
		return locRet / iteRound.getArea();
	}
	//计算适应度函数值的方法
	public void calOutsideFitness() {
		calinsideFitness();
		for(int i = 0;i < randomPopList.size();i++) {
			double e_trs = 0.0;
			double t_trs = 0.0;
			double e_inv = 0.0;
			double e_cmp = 0.0;
			double e = 0.0;
			double t = 0.0;
			double t_cmp = 0.0;
			float spt = 0.0f;
			double lbd = 0.0;
			double lbn = 0.0;
			ArrayList<Integer> temp = new ArrayList<Integer>();
			ArrayList<Double> tempE = new ArrayList<Double>();
			for(int j = 0;j < randomPopList.get(i).deviceId.size();j++) {
				for(int m = 0;m < randomPopList.get(i).deviceId.size();m++) {
					//temp.add(randomPopList.get(i).deviceId.get(m));
					if(randomPopList.get(i).deviceId.get(m).equals(randomPopList.get(i).deviceId.get(j))
							&& m != randomPopList.get(i).deviceId.size() - 1
							&& randomPopList.get(i).deviceId.get(m) != randomPopList.get(i).deviceId.get(m+1)
							&& !temp.contains(randomPopList.get(i).deviceId.get(j))) {
						e_trs += (0.2 * service.get(serviceReq.get(m)).getDt() 
								/ calRij(device.get(randomPopList.get(i).deviceId.get(m)), device.get(randomPopList.get(i).deviceId.get(m + 1))));
						t_trs += service.get(serviceReq.get(m)).getDt() 
								/ calRij(device.get(randomPopList.get(i).deviceId.get(m)), device.get(randomPopList.get(i).deviceId.get(m + 1)));
					}
					if(randomPopList.get(i).deviceId.get(m).equals(randomPopList.get(i).deviceId.get(j)) && 
							!temp.contains(randomPopList.get(i).deviceId.get(j))) {
						e_inv += service.get(serviceReq.get(m)).getInv();
						e_cmp += Math.pow(10, -26) * device.get(randomPopList.get(i).deviceId.get(m)).getF() 
								* device.get(randomPopList.get(i).deviceId.get(m)).getF()
								* service.get(serviceReq.get(m)).getCr();
						t_cmp += service.get(serviceReq.get(m)).getCr()
								/ device.get(randomPopList.get(i).deviceId.get(m)).getF();
					}
				}
				
				if(!temp.contains(randomPopList.get(i).deviceId.get(j))){
					temp.add(randomPopList.get(i).deviceId.get(j));
					e += (e_inv + e_cmp + e_trs);
					tempE.add((e_inv + e_cmp + e_trs));
					e_inv = 0.0;
					e_cmp = 0.0;
					e_trs = 0.0;
				}
			}
			for(int x = 0;x < temp.size();x++) {
				spt += evaluateLoc(device.get(temp.get(x)), iteRound);
				lbd += (device.get(temp.get(x)).getRsd() - tempE.get(x)) / device.get(temp.get(x)).getRsd();
			}
			t += t_cmp + t_trs;
			spt = spt / temp.size();
			lbn = lbd / temp.size();
			popList.get(i).fitness[0] = e;
			popList.get(i).fitness[1] = t;
			popList.get(i).fitness[2] = spt;
			popList.get(i).fitness[3] = lbn;
		}
		/*System.out.println();
		for(int i = 0;i < randomPopList.size();i++) {
			for(int j = 0;j < randomPopList.get(i).deviceId.size();j++) {
				System.out.print(randomPopList.get(i).deviceId.get(j)+" ");
			}
			System.out.println();
			for(int j = 0;j < 4;j++) {
				System.out.print(init.getPopList().get(i).fitness[j]+" ");
			}
			System.out.println();
		}*/
//		System.out.println("测试初始种群的FItness值是否正确");
//		for(int i = 0;i < randomPopList.size();i++) {
//			for(int j = 0;j < 4;j++) {
//				System.out.print(init.getPopList().get(i).fitness[j]+" ");
//			}
//			System.out.println();
//		}
//		copyPopListToPrePopList();
//		System.out.println("测试复制的种群");
//		for(int i = 0;i < randomPopList.size();i++) {
//			for(int j = 0;j < 3;j++) {
//				System.out.print(init.getPrePopList().get(i).fitness[j]+" ");
//			}
//			System.out.println();
//		}
	}
	/*
	* 将当代种群的相关数据复制到prePopList中，供下次迭代时使用
	*/
	public void copyPopListToPrePopList() {
		for (int i = 0; i < popList.size(); i++) {
			prePopList.set(i, (Chromsome)popList.get(i).clone());
		}
	}	
	/*
	 * 判断两个个体p和q之间的pareto支配关系：tag=1 表示p支配q; tag=2 表示q支配p; tag=0
	 * 表示p和q之间等价，没有支配关系.
	 */
	public int isDominate(Chromsome p, Chromsome q) { 
		int tag = 0;// tag=0表示p和q之间等价，没有支配关系
		int index1 = 0;
		for (int i = 0; i < 4; i++) { // 判断p支配q否
			if (p.fitness[i] > q.fitness[i]) {
				index1 += 1;
			}
		}
		if (index1 == 4) {
			tag = 1;// 表示p支配q
		}

		int index2 = 0;
		for (int i = 0; i < 4; i++) { // 判断q支配p否
			if (q.fitness[i] > p.fitness[i]) {
				index2 += 1;
			}
		}
		if (index2 == 4) {
			tag = 2;// 表示q支配p
		}
		return tag;
	}
	/*
	 * 对list列表中的个体按照fintness值从小到大进行排序(快速排序)
	 */
	public ArrayList<Chromsome> getSortList(ArrayList<Chromsome> list,int numObjective) {
		Chromsome temp = new Chromsome();
		for (int i = 0; i < list.size(); i++) {
			for (int j = i + 1; j < list.size(); j++) {
				if (list.get(i).fitness[numObjective] > list.get(j).fitness[numObjective]) {
					temp = list.get(i);
					list.set(i, list.get(j));
					list.set(j, temp);
				}
			}
		}
		return list;
	}
	/*
	 * 快速非支配排序
	 */
	public ArrayList<ArrayList<Chromsome>> fastNondominateSort() {
		dominatingMap = new HashMap<Chromsome,ArrayList<Chromsome>>();
		// 第一个pareto等级集合
		ArrayList<Chromsome> firstParetoRankSet = new ArrayList<Chromsome>();
		// 所有pareto等级集合
		ArrayList<ArrayList<Chromsome>> paretoRankSetList = new ArrayList<ArrayList<Chromsome>>();// pareto分级的解集
		ArrayList<Chromsome> unionPopList = new ArrayList<Chromsome>();
		for (int i = 0; i < 2 * popList.size(); i++) {
			Chromsome chromsome = new Chromsome(7);
			unionPopList.add(chromsome);
		}
		// 将prePopList中的元素复制到unionPopList中
		for (int i = 0; i < prePopList.size(); i++) {
			unionPopList.set(i, (Chromsome)prePopList.get(i).clone());
		}

		// 将PopList中的元素复制到unionPopList中
		for (int i = popSize; i < (2 * popList.size()); i++) {
			unionPopList.set(i,(Chromsome)popList.get(i-popList.size()).clone());
		}
		/*System.out.println("测试一下联合的有无适应度");
		for (int i = 0; i < 2 * prePopList.size(); i++) {
			//unionPopList.set(i, (Chromsome)init.getPrePopList().get(i).clone());
			for (int k = 0; k < 4; k++) {
				System.out.print(unionPopList.get(i).fitness[k] + " "); 
			}
			System.out.println();
		}*/
		// 每个个体和所有个体比较,判断支配
		for (int i = 0; i < unionPopList.size(); i++) {
			ArrayList<Chromsome> tempDominating = new ArrayList<Chromsome>();
			for (int j = 0; j < unionPopList.size(); j++) {
				int tag = isDominate(unionPopList.get(i), unionPopList.get(j));
				if (tag == 1) {
					tempDominating.add(unionPopList.get(j));// 记录当前个体支配的解
				} else if (tag == 2) {
					unionPopList.get(i).numDominated += 1;// 支配当前个体的个体数目加1
				}
			}
			dominatingMap.put(unionPopList.get(i), tempDominating);
			// 如果当前解的numDominated属性等于0，即属于支配最前沿，则将其加入第一个pareto等级集合，并将其paretoRank值设置为1
			if (unionPopList.get(i).numDominated == 0) { 
				unionPopList.get(i).setParetoRank(1);
				firstParetoRankSet.add(unionPopList.get(i));
			}
		}
		paretoRankSetList.add(firstParetoRankSet);
//		System.out.println("pareto front 集合的大小： " + firstParetoRankSet.size());
//		for(int i = 0 ; i < firstParetoRankSet.size(); i++) {
//			for (int j = 0 ; j < init.getNumService(); j++) {
//				System.out.print(firstParetoRankSet.get(i).deviceId.get(j) + " ");
//			}
//			System.out.println();
//		}
			
		int rank = 0;
		while (paretoRankSetList.get(rank).size() > 0
				&& paretoRankSetList.get(rank) != null) {
			// 用于储存下一前沿的pareto集合
			ArrayList<Chromsome> paretoRankSet = new ArrayList<Chromsome>();
			// 依次处理当前pareto等级里的所有个体
			for (int j = 0; j < paretoRankSetList.get(rank).size(); j++) {
				//Chromsome currentChromsome = paretoRankSetList.get(rank).get(j);
				ArrayList<Chromsome> currentDList = dominatingMap.get(paretoRankSetList.get(rank).get(j));
				if (currentDList.size() > 0 && currentDList != null) {
					for (int k = 0; k < currentDList.size(); k++) {
						Chromsome dominatedChromsome = currentDList.get(k);
						dominatedChromsome.numDominated -= 1;
						// 若当前个体属于下一个前沿，修改个体的paretoRank属性值,并将其加入相应的pareto等级集合
						if (dominatedChromsome.numDominated == 0) {
							dominatedChromsome.setParetoRank(rank + 1);
							paretoRankSet.add(dominatedChromsome);
						}
					}
				}
			}
			if (paretoRankSet != null && paretoRankSet.size() > 0) {// 如果paretoRankSet不为空，将其加入paretoRankSetList中；否则，说明所有个体都已经处理完毕，退出循环
				paretoRankSetList.add(paretoRankSet);// 将新生成的paretoRankSet加入paretoRankSetList中供后面计算crowdingDistance时使用
				rank += 1;
			} else {
				break;
			}
		}
//		for(int i = 0 ; i < paretoRankSetList.size(); i++) {
//			for (int j = 0 ; j < paretoRankSetList.get(i).size(); j++) {
//				for(int m = 0 ; m < init.getNumService() ; m++) {
//					System.out.print(paretoRankSetList.get(i).get(j).deviceId.get(m) + " ");
//				}
//				System.out.println();
//			}
//			System.out.println();
//		}
		return paretoRankSetList;
	}
	
	/*
	 *  拥挤距离计算
	 */
	public void crowdingDistance(ArrayList<ArrayList<Chromsome>> paretoRankSetList) {
		for (int i = 0; i < paretoRankSetList.size(); i++) {
			ArrayList<Chromsome> paretoRankSet = new ArrayList<Chromsome>();
			paretoRankSet = paretoRankSetList.get(i);

			if (paretoRankSet.size() > 0 && paretoRankSet != null) {
				//按照目标依次计算个体的crowdingDistance值
				for (int j = 0; j < 4; j++) {
					//对paretoRankSet中的个体按照fitness从小到大进行排序
					paretoRankSet = getSortList(paretoRankSet, j);
					if (paretoRankSet.size() == 1) {// paretoRankSet中只有一个个体
						paretoRankSet.get(0).crowdingDistance += 100.0;
					} else if (paretoRankSet.size() == 2) {// paretoRankSet中有两个个体
						paretoRankSet.get(0).crowdingDistance += 100.0;
						paretoRankSet.get(paretoRankSet.size() - 1).crowdingDistance += 1000000.0;

					} else {// paretoRankSet个体数大于2
						double minFitness = paretoRankSet.get(0).fitness[j];// 目标值的最小值
						double maxFitness = paretoRankSet.get(paretoRankSet.size() - 1).fitness[j];// 目标值的最大值

						paretoRankSet.get(0).crowdingDistance += 100.0;// 将第1个和最后1个个体的Distance设置为无穷大，以便保存处在边界上的个体，为了便于计算设置为10000.0
						paretoRankSet.get(paretoRankSet.size() - 1).crowdingDistance += 1000000.0;

						for (int k = 1; k < paretoRankSet.size() - 2; k++) {// 依次计算其他个体的crowdingDistance
							paretoRankSet.get(k).crowdingDistance += ((paretoRankSet.get(k + 1).fitness[j] - paretoRankSet.get(k - 1).fitness[j]) / (maxFitness - minFitness));
						}
					}
				}
			}
		}

		// 在一个paretoRankSet中根据个体的crowdingDistance从大到小对个体进行排序
		for (int i = 0; i < paretoRankSetList.size(); i++) {
			ArrayList<Chromsome> paretoRankSet = paretoRankSetList.get(i);
			// 对上面得到的paretoRankSet中的个体根据其crowdingDistance从大到小进行排序
			for (int j = 0; j < paretoRankSet.size(); j++) {
				Chromsome temp = new Chromsome();
				for (int k = j + 1; k < paretoRankSet.size(); k++) {
					if (paretoRankSet.get(j).crowdingDistance < paretoRankSet.get(k).crowdingDistance) {
						temp = paretoRankSet.get(j);
						paretoRankSet.set(j, paretoRankSet.get(k));
						paretoRankSet.set(k, temp);
					}
				}
			}
		}
	}
	
	public void makeNewPopulation() {
		// 对两代种群快速非支配排序
		ArrayList<ArrayList<Chromsome>> paretoRankSetListAll = fastNondominateSort();
		// 对每个pareto前沿计算拥挤距离并排序
		crowdingDistance(paretoRankSetListAll);
		// 记录加入到popList中的个体总数
		int numChromsome = 0;
		// popList的下标
		int popListIndex = 0;
		//取出前面n个个体作为下一代种群
			for (int i = 0; i < paretoRankSetListAll.size(); i++) {
				ArrayList<Chromsome> paretoRankSet = paretoRankSetListAll.get(i);
				numChromsome += paretoRankSet.size();

				if (numChromsome < popSize && popListIndex < popSize) {
					for (int m = 0; m < paretoRankSet.size(); m++) {
						popList.set(popListIndex++, paretoRankSet.get(m));
						// popList.get(popListIndex++).setBinChrom(paretoRankSet.get(j).getBinChrom());
					}
				} else {
					int overNum = numChromsome - popSize;
					for (int m = 0; m < paretoRankSet.size() - overNum && popListIndex < popSize; m++) {
						popList.set(popListIndex++, paretoRankSet.get(m));
						 
						// popList.get(popListIndex++).setBinChrom(paretoRankSet.get(j).getBinChrom());
					}
				}
			}
		
		/*System.out.println("新一代");
		for(int i = 0;i < popList.size();i++) {
				for(int j = 0;j < numService;j++) {
					System.out.print(popList.get(i).deviceId.get(j)+" ");
				}
				System.out.println();
				for(int x = 0;x < 4;x++) {
					System.out.print(popList.get(i).fitness[x]+" ");
			}
			System.out.println();
			System.out.println(0.4*popList.get(i).fitness[0] + 0.4*popList.get(i).fitness[1] + popList.get(i).fitness[2]*0.1+popList.get(i).fitness[3]*0.1);
			System.out.println();
		}*/
	}
	
	/*
	 * 1 对种群中的个体进行随机配对并将配对结果保存到index中 2 对配对后的个体执行单点交叉，按照index中的次序相邻两个个体进行交叉 3
	 * 交叉后的个体直接代替父代个体进入种群 4 未对最优个体采取任何保护策略
	 */
	public void crossoverOperator() {
		Random random = new Random();
		int point;// 配对过程的中间变量
		int temp;// 交换过程中的辅助变量

		// 对种群中的个体进行随机配对
		int[] index = new int[popList.size()];
		for(int i = 0; i < popList.size() ; i++) {
			index[i] = i;
			//System.out.println(index[i] + " ");
		}
		for (int i = 0; i < popList.size(); i++) {
			point = random.nextInt(popList.size());
			temp = index[i];
			index[i] = index[point];
			index[point] = temp;
		}
		for (int i = 0; i < popList.size(); i += 2) {
			double pro = random.nextDouble();
			if (pro < pc){
				point = random.nextInt(numService);// 交叉点
				// 交换index[i]和index[i+1]在交叉点后面的个体
				for (int j = point; j < numService; j++) {
					temp =popList.get(index[i]).deviceId.get(j);
					popList.get(index[i]).deviceId.set(j,popList.get(index[i+1]).deviceId.get(j));
					popList.get(index[i+1]).deviceId.set(j, temp);
				}
			}
		}
	}
	// 按概率pm进行单点变异，随机选择一个变异点，然后在同类服务集合里面随机找一个设备替代它
	public void mutationOperator() {
		Random random = new Random();
		int deviceId1 = 0;
		int deviceId2 = 0;
		int serviceId1 = 0;
		double pro = 0.0;// 生成一个随机数用于根据变异概率pm判断当前个体是否参与变异
		int index1 = 0; //生成一个随机数用于染色体上的某一个点发生变异
		int index2 = 0; //生成一个随机数，在服务集合里面随机找一个不同的设备，使得该设备发生变化
		for (int i = 0; i < popList.size(); i++) {//遍历染色体
			pro = random.nextDouble();
			if (pro < pm) {
				index1 = random.nextInt(numService);//产生的随机数的范围是：服务数量-1
				deviceId1 = popList.get(i).deviceId.get(index1);//记录下来此时的设备Id，从而获得该设备ID的对应服务。
				serviceId1 = serviceReq.get(index1);//记录下来此时的服务Id，从该服务类的集合中找设备
				//从该服务Id中随机产生一个设备号
				//若随机产生的设备号和该设备号不同时，替换他，否则，重新产生随机数。
				index2 = random.nextInt(service.get(serviceId1).getDevice().size());
				deviceId2 = service.get(serviceId1).getDevice().get(index2).getId();
				while(deviceId1 == deviceId2) {
					index2 = random.nextInt(service.get(serviceId1).getDevice().size());
					deviceId2 = service.get(serviceId1).getDevice().get(index2).getId();
				}
				popList.get(i).deviceId.set(index1, deviceId2);			
			}
		}
//		System.out.println("测试开始");
//		for(int i = 0; i < init.getPopList().size(); i++ ) {
//			for (int j = 0; j < init.getNumService(); j++) {
//				System.out.print(init.getPopList().get(i).deviceId.get(j) + " ");
//			}
//			System.out.println();
//		}
	} 	
	/*
	 * 运行算法逻辑,图一二三
	 * */
	////////////////////////////////////////
	/*public ArrayList<Double> runNSGA2(){
		ArrayList<Double> U_max2 = new ArrayList<>();
		//计算当前种群的适应度
		calOutsideFitness();
		// 将种群保存起来供下次迭代
		copyPopListToPrePopList();
		int gen = 0;
		while (gen < maxGen) {
			//进入下一代
			gen += 1;
			//交叉变异操作
			crossoverOperator();
			mutationOperator();
			// 计算目前适应度
			calOutsideFitness();
			//找出支配高，拥挤距离大的个体形成新的种群
			//System.out.println("第"+gen+"代");
			makeNewPopulation();
			//前三幅图需要用到的循环
			if(gen >= 10) {
				System.out.println("第"+gen+"代随机迁移外部最优");
				double max_u= 0.0;
				max_u = chooseOutside();
				U_max2.add(max_u);
			}
			//前三幅图需要用到的循环
			//将选取出来的种群保存到上一代中
			copyPopListToPrePopList();
		}		
		return U_max2;
	}*/
	/*
	 * 图四和图五、图八和图九需要用到的算法运行逻辑
	 */
	public Map runNSGA2(){
		Map randomResult = new HashMap();
		double max_u= 0.0;
		//计算当前种群的适应度
		calOutsideFitness();
		// 将种群保存起来供下次迭代
		copyPopListToPrePopList();
		int gen = 0;
		while (gen < maxGen) {
			//进入下一代
			gen += 1;
			//System.out.println(gen);
			//交叉变异操作
			crossoverOperator();
			mutationOperator();
			// 计算目前适应度
			calOutsideFitness();
			//找出支配高，拥挤距离大的个体形成新的种群
			//System.out.println("第"+gen+"代");
			makeNewPopulation();
			//将选取出来的种群保存到上一代中
			copyPopListToPrePopList();
		}
		randomResult = chooseOutside();
		return randomResult;
	}

	public Map chooseOutside() {
		ArrayList<Chromsome> bestSolution = new ArrayList<Chromsome>();
		double min_e = popList.get(0).fitness[0];
		double max_e = popList.get(0).fitness[0];
		double min_t = popList.get(0).fitness[1];
		double max_t = popList.get(0).fitness[1];
		double min_spt = popList.get(0).fitness[2];
		double max_spt = popList.get(0).fitness[2];
		double min_lbn = popList.get(0).fitness[3];
		double max_lbn = popList.get(0).fitness[3];
		double[] u_e = new double[popList.size()];
		double[] u_t = new double[popList.size()];
		double[] u_spt = new double[popList.size()];
		double[] u_lbn = new double[popList.size()];
		double[] u = new double[popList.size()];
		for(int i = 1;i < popList.size();i++) {
			if(popList.get(i).fitness[0] < min_e) {
				min_e = popList.get(i).fitness[0];
			}
			if(popList.get(i).fitness[0] > max_e) {
				max_e = popList.get(i).fitness[0];
			}
			if(popList.get(i).fitness[1] < min_t) {
				min_t = popList.get(i).fitness[1];
			}
			if(popList.get(i).fitness[1] > max_t) {
				max_t = popList.get(i).fitness[1];
			}
			if(popList.get(i).fitness[2] < min_spt) {
				min_spt = popList.get(i).fitness[2];
			}
			if(popList.get(i).fitness[2] > max_spt) {
				max_spt = popList.get(i).fitness[2];
			}
			if(popList.get(i).fitness[3] < min_lbn) {
				min_lbn = popList.get(i).fitness[3];
			}
			if(popList.get(i).fitness[3] > max_lbn) {
				max_lbn = popList.get(i).fitness[3];
			}
		}
		for(int i = 0;i < popList.size();i++) {
			if(min_e == max_e) {
				u_e[i] = 1;
			}
			else {
				u_e[i] = (max_e - popList.get(i).fitness[0]) / (max_e - min_e);
			}
			if(min_t == max_t) {
				u_t[i] = 1;
			}
			else {
				u_t[i] = (max_t - popList.get(i).fitness[1]) / (max_t - min_t);
			}
			if(min_spt == max_spt) {
				u_spt[i] = 1;
			}
			else {
				u_spt[i] = (max_spt - popList.get(i).fitness[2]) / (max_spt - min_spt);
			}
			if(min_lbn == max_lbn) {
				u_lbn[i] = 1;
			}
			else {
				u_lbn[i] = (max_lbn - popList.get(i).fitness[3]) / (max_lbn - min_lbn);
			}	
			u[i] = 0.4 * u_e[i] + 0.4 * u_t[i] - 0.1 * u_spt[i] - 0.1 * u_lbn[i]; 
		}
		double max_u = u[0];
		int point = 0;
		for(int i = 1;i < popList.size();i++) {
			if(u[i] > max_u) {
				point = i;
				max_u = u[i];
 			}
		}
		bestSolution.add(popList.get(point));
		//System.out.println("外部最优");
//		for(int i = 0;i < bestSolution.size();i++) {
//			for(int j = 0;j < bestSolution.get(i).deviceId.size();j++) {
//				System.out.print(bestSolution.get(i).deviceId.get(j)+" ");
//			}
//			System.out.println();
//		}
//		System.out.println("对应策略：");
//		for(int i = 0;i < bestSolution.size();i++) {
//			for(int j = 0;j < bestSolution.get(i).mgDeviceId.size();j++) {
//				System.out.print(bestSolution.get(i).mgDeviceId.get(j)+" ");
//			}
//			System.out.println();
//		}
		//图二需要计算的每一代的e_min
		//图三需要计算的每一代的能量的方差e_variance
		double e_variance=0.0;//能量的方差
		double e_min = 50.0;//能量最小值
		double e_eve = 50.0;// 每个设备的能量值
		double temp = 0.0;// 中间变量,记录迁移能耗
		double e_trs = 0.0;// 服务组合某个设备的传输能耗
		double e_inv = 0.0;// 服务组合某个设备的调用能耗
		double e_cmp = 0.0;// 服务组合某个设备的计算能耗
		ArrayList<Double> E_eve = new ArrayList<>();// 记录当前最优染色体的每个设备的能量
		// 计算最优个体迁移所消耗的能量和服务组合所消耗的能量。
		// 首先计算该服务组合的迁移能耗
		for (int i = 0; i < bestSolution.size(); i++) {
			for (int j = 0; j < bestSolution.get(i).deviceId.size(); j++) {
				for (int m = 0; m < bestSolution.get(i).deviceId.size(); m++) {
					if (bestSolution.get(i).deviceId.get(j).equals(bestSolution.get(i).mgDeviceId.get(m)) && j != m) {
						temp += 0.2 * service.get(serviceReq.get(m)).getWkd() * 1048576
								/ calRij(device.get(bestSolution.get(i).deviceId.get(j)),
										device.get(bestSolution.get(i).deviceId.get(m)));
					}
					if (bestSolution.get(i).deviceId.get(j).equals(bestSolution.get(i).mgDeviceId.get(m))
							&& m != bestSolution.get(i).deviceId.size() - 1
							&& bestSolution.get(i).mgDeviceId.get(m) != bestSolution.get(i).mgDeviceId.get(m + 1)) {
						e_trs += (0.2 * service.get(serviceReq.get(m)).getDt()
								/ calRij(device.get(bestSolution.get(i).mgDeviceId.get(m)),
										device.get(bestSolution.get(i).mgDeviceId.get(m + 1))));
					}
					if (bestSolution.get(i).deviceId.get(j).equals(bestSolution.get(i).mgDeviceId.get(m))) {
						e_inv += service.get(serviceReq.get(m)).getInv();
						e_cmp += Math.pow(10, -26) * device.get(bestSolution.get(i).deviceId.get(j)).getF()
								* device.get(bestSolution.get(i).deviceId.get(j)).getF()
								* service.get(serviceReq.get(m)).getCr();
					}
				}
				e_eve = e_eve - temp - e_trs - e_inv - e_cmp;
				E_eve.add(e_eve);
				e_eve = 50.0;
				e_trs = 0.0;
				e_inv = 0.0;
				e_cmp = 0.0;
				temp = 0.0;
			}
		}
		e_variance = getVariance(E_eve);
//		System.out.println("测试一下E_eve的值是否正确");
//		for (int i = 0; i < E_eve.size(); i++) {
//				System.out.print(E_eve.get(i) + " ");
//		}
//		System.out.println();
		for (int i = 0; i < E_eve.size(); i++) {
			if(E_eve.get(i) < e_min) {
				e_min = E_eve.get(i);
			}
		}
		//画图一时，需要返回的值是max_u
		//return max_u;
		//画图二时，需要返回的值是bestSolution中的剩余能量最小的设备的能量值e_min。
		//return e_min;
		//画图三时，需要返回的值是bestSolution的能量方差
		//return e_variance;
		//画图四和图五、图八和图九时，需要分别输出服务组合能耗conErg和服务组合时间conTime
		double conErg = 0.0;
		double conTime = 0.0;
		//服务迁移时间和用户数量之间的关系，这里定义一个内层的消耗时间变量
		double insideconTime = 0.0;
		double insideconErg = 0.0;
		//服务组合消耗的总能量应是迁移能耗+服务组合能耗，即内层的fitness值和外层的fitness值之和
		//System.out.println("输出内层的适应度函数值： " + bestSolution.get(0).insidefitness[1]);
		//System.out.println("输出外层的适应度函数值： " + bestSolution.get(0).fitness[0]);
		for(int i = 0; i < bestSolution.size(); i++) {
			conErg = bestSolution.get(i).fitness[0];
			conTime = bestSolution.get(i).fitness[1];
			insideconTime = bestSolution.get(i).insidefitness[0];
			insideconErg = bestSolution.get(i).insidefitness[1];
		}
		//System.out.println("conErgde的值是" + conErg);
		//return conErg;
		//return conTime;
		//return insideconTime;
		Map randomResult = new HashMap();

		randomResult.put("bestSolution",bestSolution);
		randomResult.put("conErg",conErg);
		randomResult.put("conTime",conTime);
		return randomResult;
		
	}

	public double getVariance(ArrayList<Double> x) {
		int m = x.size();
		double sum = 0;
		for (int i = 0; i < m; i++) {// 求和
			sum += x.get(i);
		}
		double dAve = sum / m;// 求平均值
		double dVar = 0;
		for (int i = 0; i < m; i++) {// 求方差
			dVar += (x.get(i) - dAve) * (x.get(i) - dAve);
		}
		return dVar / m;
	}

	// 主函数，测试算法输出
	/*public static void main(String[] args) {
		ArrayList<Chromsome> List= new ArrayList<>();
		int maxGen = 200; //最大进化代数
		double pc = 0.6;// 交叉概率
		double pm = 0.6;// 变异概率
		RandomMig R1 = new RandomMig(maxGen,pc, pm);
//		InitPopulation T1 = R1.getInit();
//		System.out.println("初始化产生的染色体长度是：" + T1.getNumService());// 染色体长度
		//R1.calOutsideFitness();
		//R1.fastNondominateSort();
		//R1.makeNewPopulation();
		List = R1.runNSGA2();
		System.out.println();
		R1.selectPopList(List);
	}*/

}