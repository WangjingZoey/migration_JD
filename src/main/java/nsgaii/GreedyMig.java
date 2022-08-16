package nsgaii;
/*贪婪算法
 * 三个约束条件
 * 贪婪策略：max（0.5*剩余能量+0.5*剩余能量）。
 * 
 */
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import networkInit.CreateNetwork;
import networkInit.DeviceInit;
import networkInit.Round;
import variable.Variable;

public class GreedyMig {
	InitPopulation init;

	private int maxGen; // 最大进化代数
	// private int numService;
	private int gen;// 当前进化代数
	private int lowBand;
	private int highBand;
	double pc;// 交叉概率

	double pm;// 变异概率
	Random random = new Random();
	ArrayList<Chromsome> greedyPopList;
	private float[] Cache;// 存储每个设备剩余容量数组
	private float[] Energy;// 存储每个设备剩余能量数组
	private float[] Fitness;// 存储每个设备的最大剩余能量和容量之和
	private int[] Migration;// 存储设备号列表进行复制操作
	HashMap<Chromsome, ArrayList<Chromsome>> dominatingMap;
	CreateNetwork network;

	public GreedyMig(InitPopulation init, int maxGen, double pc, double pm, CreateNetwork network,int lowBand,int highBand) {
		this.init = init;
		this.maxGen = maxGen;
		this.pc = pc;
		this.pm = pm;
		this.network = network;
		// this.numService = numService;
		Cache = new float[init.getNumService()];
		Energy = new float[init.getNumService()];
		Fitness = new float[init.getNumService()];
		Migration = new int[init.getNumService()];
		this.lowBand = lowBand;
		this.highBand = highBand;
	}

	// 约束条件一，每个设备上的服务数量不超过五个
	public boolean JdgServiceNum(ArrayList<Integer> x) {
		int a = 1;
		for (int i = 0; i < x.size() - 1; i++) {// for
			// System.out.print(x.get(i) + " ");//get():获取指定索引处的值'
			for (int j = i + 1; j < x.size(); j++) {
				if (x.get(i) == x.get(j))
					a++;
			}
			if (a > 5)
				return true;
			a = 1;
		}
		return false;
	}

	// 复制染色体集合
	public void copyPopListToGreedy() {
		// init.initPopulation();// 生成初始种群
		greedyPopList = new ArrayList<Chromsome>();
		// 复制PopList；
		for (int i = 0; i < init.getPopList().size(); i++) {
			Chromsome chromsome = new Chromsome(7);
			greedyPopList.add(chromsome);
			greedyPopList.set(i, (Chromsome)init.getPopList().get(i).clone());
		}
	}

	// 初始化Cache，Energy数组
	public void initArray() {
		for (int i = 0; i < init.getPopList().size(); i++) {
			for (int j = 0; j < init.getNumService(); j++) {
				Cache[j] = 8.0f - init.getService().get(init.getServiceReq().get(j)).getWkd();
				Energy[j] = init.getDevice().get(init.getPopList().get(i).deviceId.get(j)).getRsd();
			}
		}
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
		return rij;
	}
	// 贪婪迁移算法实现逻辑
	public void greedyMigration() {
		copyPopListToGreedy();
		float maxNum = 0.0f;
		int IdDevice = 0;
		float e_trs = 0.0f;
		int a = 0;
		boolean count;
		// System.out.println();
//		for (int i = 0; i < init.getPopList().size(); i++) {
//			for (int j = 0; j < init.getNumService(); j++) {
//				Migration[j] = init.getPopList().get(i).deviceId.get(j);
//				greedyPopList.get(i).deviceId.add(j, Migration[j]);
//			}
//		}
//		System.out.println("测试开始");
//		for (int i = 0; i < init.getPopList().size(); i++) {
//			for (int j = 0; j < init.getNumService(); j++) {
//				System.out.print(init.getPopList().get(i).deviceId.get(j) + " ");
//			}
//			System.out.println();
//		}
//		System.out.println("greedyPopList");
//		for (int i = 0; i < init.getPopList().size(); i++) {
//			for (int j = 0; j < init.getNumService(); j++) {
//				System.out.print(greedyPopList.get(i).deviceId.get(j) + " ");
//			}
//			System.out.println();
//		}
		initArray();
		for (int i = 0; i < init.getPopList().size(); i++) {
			for (int m = 0; m < init.getNumService(); m++) {// 从第一个设备开始迁移
				for (int j = 0; j < init.getNumService(); j++) {// 从第一个设备开始计算
					Fitness[j] = 0.5f * Cache[j] + 0.5f * Energy[j];// 贪婪策略，计算适应度函数值的最大值
					// System.out.println(init.getPopList().get(i).deviceId.get(j) + "号设备的容量是：" +
					// Cache[j]);
					// System.out.println(init.getPopList().get(i).deviceId.get(j) + "号设备的能量是："
					// +Energy[j]);
					// System.out.print(init.getPopList().get(i).deviceId.get(j) + "号设备的贪婪函数值是："
					// +Fitness[j]+ " ");
				}
				// System.out.println();
				// 计算适应度函数值里的最大值
				for (int n = 0; n < Fitness.length; n++) {
					if (Fitness[n] > maxNum) {
						maxNum = Fitness[n];
						IdDevice = n;
					}
				}
				maxNum = 0.0f;
				// 将第m处的设备迁移到第IdDevice处设备上，
				a = init.getPopList().get(i).deviceId.get(IdDevice);
				greedyPopList.get(i).deviceId.set(m, a);
				// greedyPopList.get(i).deviceId.set(m,init.getPopList().get(i).deviceId.get(IdDevice));
				// 更新每个设备的剩余容量和剩余能量
				Cache[IdDevice] = Cache[IdDevice] - init.getService().get(init.getServiceReq().get(m)).getWkd();
				Cache[m] = Cache[m] + init.getService().get(init.getServiceReq().get(m)).getWkd();
				e_trs = (float) (0.2 * init.getService().get(init.getServiceReq().get(m)).getWkd() * 1048576
						/ calRij(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)),
								init.getDevice().get(init.getPopList().get(i).deviceId.get(IdDevice))));
				Energy[m] = Energy[m] - e_trs;
				count = JdgServiceNum(greedyPopList.get(i).deviceId);
				// init.getService().get(init.getServiceReq().get(m)).getWkd().set()
				while (Cache[IdDevice] < 0 || Energy[m] < 0 || count == true) {
					initArray();
					for (int p = 0; p < init.getNumService(); p++) {// 从第一个设备开始迁移
						for (int q = 0; q < init.getNumService(); q++) {// 从第一个设备开始计算
							Fitness[q] = 0.5f * Cache[q] + 0.5f * Energy[q];// 贪婪策略，计算适应度函数值的最大值
						}
						for (int n = 0; n < Fitness.length; n++) {
							if (Fitness[n] > maxNum) {
								maxNum = Fitness[n];
								IdDevice = n;
							}
						}
						maxNum = 0.0f;
						// 将第p处的设备迁移到第IdDevice号设备上
						a = init.getPopList().get(i).deviceId.get(IdDevice);
						greedyPopList.get(i).deviceId.set(p, a);
						// greedyPopList.get(i).deviceId.set(m,init.getPopList().get(i).deviceId.get(IdDevice));
						// 更新每个设备的剩余容量和剩余能量
						Cache[IdDevice] = Cache[IdDevice] - init.getService().get(init.getServiceReq().get(p)).getWkd();
						Cache[p] = Cache[p] + init.getService().get(init.getServiceReq().get(p)).getWkd();
						e_trs = (float) (0.2 * init.getService().get(init.getServiceReq().get(p)).getWkd() * 1048576
								/ calRij(init.getDevice().get(init.getPopList().get(i).deviceId.get(p)),
										init.getDevice().get(init.getPopList().get(i).deviceId.get(IdDevice))));
						Energy[p] = Energy[p] - e_trs;
						count = JdgServiceNum(greedyPopList.get(i).deviceId);
					}
				}
			}

		}
		for(int i = 0; i < greedyPopList.size(); i++) {
			for (int j=0 ; j < greedyPopList.get(i).deviceId.size(); j++) {
				init.getPopList().get(i).mgDeviceId.set(j,greedyPopList.get(i).deviceId.get(j));
		    }
		}
	}

	// 计算迁移时间和迁移的能耗
	public void calinsideFitness() {
		greedyMigration();
		double t = 0.0;
		double e = 0.0;
		double cst = 0.0;
		double cbf = 0.0;
		for (int i = 0; i < init.getPopList().size(); i++) {
			double[] fitness = new double[3];
			for (int m = 0; m < init.getNumService(); m++) {
				for (int n = 0; n < init.getNumService(); n++) {
					if (init.getPopList().get(i).deviceId.get(m).equals(greedyPopList.get(i).deviceId.get(n))) {
						cst += init.getService().get(init.getServiceReq().get(n)).getWkd();
					}
				}
				cbf += (init.getDevice().get(init.getPopList().get(i).deviceId.get(m)).getRsd() - cst)
						/ init.getDevice().get(init.getPopList().get(i).deviceId.get(m)).getRsd();
				cst = 0.0;
				if (init.getPopList().get(i).deviceId.get(m) != greedyPopList.get(i).deviceId.get(m)) {
					t += init.getService().get(init.getServiceReq().get(m)).getWkd() * 1048576
							/ calRij(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)),
									init.getDevice().get(greedyPopList.get(i).deviceId.get(m)))
							+ init.getService().get(init.getServiceReq().get(m)).getCr()
									/ init.getDevice().get(greedyPopList.get(i).deviceId.get(m)).getF();
					e += 0.2 * init.getService().get(init.getServiceReq().get(m)).getWkd() * 1048576
							/ calRij(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)),
									init.getDevice().get(greedyPopList.get(i).deviceId.get(m)))
							+ 0.002 * init.getService().get(init.getServiceReq().get(m)).getCr()
									/ init.getDevice().get(greedyPopList.get(i).deviceId.get(m)).getF();
				} else {
					if (m < (init.getNumService() - 1)) {
						t += init.getService().get(init.getServiceReq().get(m)).getCr()
								/ init.getDevice().get(greedyPopList.get(i).deviceId.get(m)).getF()
								+ init.getService().get(init.getServiceReq().get(m)).getDt()
										/ calRij(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)),
												init.getDevice().get(init.getPopList().get(i).deviceId.get(m + 1)));
						e += 0.2 * init.getService().get(init.getServiceReq().get(m)).getDt()
								/ calRij(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)),
										init.getDevice().get(init.getPopList().get(i).deviceId.get(m + 1)))
								+ 0.5 * init.getService().get(init.getServiceReq().get(m)).getCr()
										/ init.getDevice().get(greedyPopList.get(i).deviceId.get(m)).getF();
					} else if (m == (init.getNumService() - 1)) {
						t += init.getService().get(init.getServiceReq().get(m)).getCr()
								/ init.getDevice().get(greedyPopList.get(i).deviceId.get(m)).getF();
						e += 0.5 * init.getService().get(init.getServiceReq().get(m)).getCr()
								/ init.getDevice().get(greedyPopList.get(i).deviceId.get(m)).getF();
					}

				}
			}
			cbf = cbf / init.getNumService();
			init.getPopList().get(i).insidefitness[0] = t;
			init.getPopList().get(i).insidefitness[1] = e;
			init.getPopList().get(i).insidefitness[2] = cbf;
			/*
			 * System.out.print("t=" + fitness[0] + " " + "e=" + fitness[1] + " " + "cbf=" +
			 * fitness[2]); System.out.println();
			 */
			t = 0.0;
			e = 0.0;
			cbf = 0.0;
		}
	}

	/*
	 * 以下是跑外层Nsga算法的流程
	 * 
	 * 
	 * 以下是实验的结果
	 */
	// 计算两个设备之间的空间距离
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

	// 计算适应度函数值的方法
	public void calOutsideFitness() {
		calinsideFitness();
		for (int i = 0; i < greedyPopList.size(); i++) {
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
			for (int j = 0; j < greedyPopList.get(i).deviceId.size(); j++) {
				for (int m = 0; m < greedyPopList.get(i).deviceId.size(); m++) {
					// temp.add(greedyPopList.get(i).deviceId.get(m));
					if (greedyPopList.get(i).deviceId.get(m).equals(greedyPopList.get(i).deviceId.get(j))
							&& m != greedyPopList.get(i).deviceId.size() - 1
							&& greedyPopList.get(i).deviceId.get(m) != greedyPopList.get(i).deviceId.get(m + 1)
							&& !temp.contains(greedyPopList.get(i).deviceId.get(j))) {
						e_trs += (0.2 * init.getService().get(init.getServiceReq().get(m)).getDt()
								/ calRij(init.getDevice().get(greedyPopList.get(i).deviceId.get(m)),
										init.getDevice().get(greedyPopList.get(i).deviceId.get(m + 1))));
						t_trs += init.getService().get(init.getServiceReq().get(m)).getDt()
								/ calRij(init.getDevice().get(greedyPopList.get(i).deviceId.get(m)),
										init.getDevice().get(greedyPopList.get(i).deviceId.get(m + 1)));
					}
					if (greedyPopList.get(i).deviceId.get(m).equals(greedyPopList.get(i).deviceId.get(j))
							&& !temp.contains(greedyPopList.get(i).deviceId.get(j))) {
						e_inv += init.getService().get(init.getServiceReq().get(m)).getInv();
						e_cmp += Math.pow(10, -26) * init.getDevice().get(greedyPopList.get(i).deviceId.get(m)).getF()
								* init.getDevice().get(greedyPopList.get(i).deviceId.get(m)).getF()
								* init.getService().get(init.getServiceReq().get(m)).getCr();
						t_cmp += init.getService().get(init.getServiceReq().get(m)).getCr()
								/ init.getDevice().get(greedyPopList.get(i).deviceId.get(m)).getF();
					}
				}

				if (!temp.contains(greedyPopList.get(i).deviceId.get(j))) {
					temp.add(greedyPopList.get(i).deviceId.get(j));
					e += (e_inv + e_cmp + e_trs);
					tempE.add((e_inv + e_cmp + e_trs));
					e_inv = 0.0;
					e_cmp = 0.0;
					e_trs = 0.0;
				}
			}
			for (int x = 0; x < temp.size(); x++) {
				spt += evaluateLoc(init.getDevice().get(temp.get(x)), network.getIteRound());
				lbd += (init.getDevice().get(temp.get(x)).getRsd() - tempE.get(x))
						/ init.getDevice().get(temp.get(x)).getRsd();
			}
			t += t_cmp + t_trs;
			spt = spt / temp.size();
			lbn = lbd / temp.size();
			init.getPopList().get(i).fitness[0] = e;
			init.getPopList().get(i).fitness[1] = t;
			init.getPopList().get(i).fitness[2] = spt;
			init.getPopList().get(i).fitness[3] = lbn;
		}
		/*
		 * System.out.println(); for(int i = 0;i < greedyPopList.size();i++) { for(int j
		 * = 0;j < greedyPopList.get(i).deviceId.size();j++) {
		 * System.out.print(greedyPopList.get(i).deviceId.get(j)+" "); }
		 * System.out.println(); for(int j = 0;j < 4;j++) {
		 * System.out.print(init.getPopList().get(i).fitness[j]+" "); }
		 * System.out.println(); }
		 */
//			System.out.println("测试初始种群的FItness值是否正确");
//			for(int i = 0;i < greedyPopList.size();i++) {
//				for(int j = 0;j < 4;j++) {
//					System.out.print(init.getPopList().get(i).fitness[j]+" ");
//				}
//				System.out.println();
//			}
//			copyPopListToPrePopList();
//			System.out.println("测试复制的种群");
//			for(int i = 0;i < greedyPopList.size();i++) {
//				for(int j = 0;j < 3;j++) {
//					System.out.print(init.getPrePopList().get(i).fitness[j]+" ");
//				}
//				System.out.println();
//			}
	}

	/*
	 * 将当代种群的相关数据复制到prePopList中，供下次迭代时使用
	 */
	public void copyPopListToPrePopList() {
		for (int i = 0; i < init.getPopList().size(); i++) {
			init.getPrePopList().set(i, (Chromsome) init.getPopList().get(i).clone());
		}
	}

	/*
	 * 判断两个个体p和q之间的pareto支配关系：tag=1 表示p支配q; tag=2 表示q支配p; tag=0 表示p和q之间等价，没有支配关系.
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
	public ArrayList<Chromsome> getSortList(ArrayList<Chromsome> list, int numObjective) {
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
	 * 
	 */
	public ArrayList<ArrayList<Chromsome>> fastNondominateSort() {
		dominatingMap = new HashMap<Chromsome, ArrayList<Chromsome>>();
		// 第一个pareto等级集合
		ArrayList<Chromsome> firstParetoRankSet = new ArrayList<Chromsome>();
		// 所有pareto等级集合
		ArrayList<ArrayList<Chromsome>> paretoRankSetList = new ArrayList<ArrayList<Chromsome>>();// pareto分级的解集
		ArrayList<Chromsome> unionPopList = new ArrayList<Chromsome>();
		for (int i = 0; i < 2 * init.getPopList().size(); i++) {
			Chromsome chromsome = new Chromsome(7);
			unionPopList.add(chromsome);
		}
		// 将prePopList中的元素复制到unionPopList中
		for (int i = 0; i < init.getPrePopList().size(); i++) {
			unionPopList.set(i, (Chromsome) init.getPrePopList().get(i).clone());
		}

		// 将PopList中的元素复制到unionPopList中
		for (int i = init.getPrePopList().size(); i < (2 * init.getPopList().size()); i++) {
			unionPopList.set(i, (Chromsome) init.getPopList().get(i - init.getPopList().size()).clone());
		}
		/*
		 * System.out.println("测试一下联合的有无适应度"); for (int i = 0; i < 2 *
		 * init.getPrePopList().size(); i++) { //unionPopList.set(i,
		 * (Chromsome)init.getPrePopList().get(i).clone()); for (int k = 0; k < 4; k++)
		 * { System.out.print(unionPopList.get(i).fitness[k] + " "); }
		 * System.out.println(); }
		 */
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
//			System.out.println("pareto front 集合的大小： " + firstParetoRankSet.size());
//			for(int i = 0 ; i < firstParetoRankSet.size(); i++) {
//				for (int j = 0 ; j < init.getNumService(); j++) {
//					System.out.print(firstParetoRankSet.get(i).deviceId.get(j) + " ");
//				}
//				System.out.println();
//			}

		int rank = 0;
		while (paretoRankSetList.get(rank).size() > 0 && paretoRankSetList.get(rank) != null) {
			// 用于储存下一前沿的pareto集合
			ArrayList<Chromsome> paretoRankSet = new ArrayList<Chromsome>();
			// 依次处理当前pareto等级里的所有个体
			for (int j = 0; j < paretoRankSetList.get(rank).size(); j++) {
				// Chromsome currentChromsome = paretoRankSetList.get(rank).get(j);
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
//			for(int i = 0 ; i < paretoRankSetList.size(); i++) {
//				for (int j = 0 ; j < paretoRankSetList.get(i).size(); j++) {
//					for(int m = 0 ; m < init.getNumService() ; m++) {
//						System.out.print(paretoRankSetList.get(i).get(j).deviceId.get(m) + " ");
//					}
//					System.out.println();
//				}
//				System.out.println();
//			}
		return paretoRankSetList;
	}

	/*
	 * 拥挤距离计算
	 * 
	 */
	public void crowdingDistance(ArrayList<ArrayList<Chromsome>> paretoRankSetList) {
		for (int i = 0; i < paretoRankSetList.size(); i++) {
			ArrayList<Chromsome> paretoRankSet = new ArrayList<Chromsome>();
			paretoRankSet = paretoRankSetList.get(i);

			if (paretoRankSet.size() > 0 && paretoRankSet != null) {
				// 按照目标依次计算个体的crowdingDistance值
				for (int j = 0; j < 4; j++) {
					// 对paretoRankSet中的个体按照fitness从小到大进行排序
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
							paretoRankSet.get(k).crowdingDistance += ((paretoRankSet.get(k + 1).fitness[j]
									- paretoRankSet.get(k - 1).fitness[j]) / (maxFitness - minFitness));
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
		// 取出前面n个个体作为下一代种群
		for (int i = 0; i < paretoRankSetListAll.size(); i++) {
			ArrayList<Chromsome> paretoRankSet = paretoRankSetListAll.get(i);
			numChromsome += paretoRankSet.size();

			if (numChromsome < init.getPrePopList().size() && popListIndex < init.getPrePopList().size()) {
				for (int m = 0; m < paretoRankSet.size(); m++) {
					init.getPopList().set(popListIndex++, paretoRankSet.get(m));
					// popList.get(popListIndex++).setBinChrom(paretoRankSet.get(j).getBinChrom());
				}
			} else {
				int overNum = numChromsome - init.getPrePopList().size();
				for (int m = 0; m < paretoRankSet.size() - overNum && popListIndex < init.getPrePopList().size(); m++) {
					init.getPopList().set(popListIndex++, paretoRankSet.get(m));

					// popList.get(popListIndex++).setBinChrom(paretoRankSet.get(j).getBinChrom());
				}
			}
		}
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
		int[] index = new int[init.getPopList().size()];
		for (int i = 0; i < init.getPopList().size(); i++) {
			index[i] = i;
			// System.out.println(index[i] + " ");
		}
		for (int i = 0; i < init.getPopList().size(); i++) {
			point = random.nextInt(init.getPopList().size());
			temp = index[i];
			index[i] = index[point];
			index[point] = temp;
		}
		for (int i = 0; i < init.getPopList().size(); i += 2) {
			double pro = random.nextDouble();
			if (pro < 0.6) {
				point = random.nextInt(init.getNumService());// 交叉点
				// 交换index[i]和index[i+1]在交叉点后面的个体
				for (int j = point; j < init.getNumService(); j++) {
					temp = init.getPopList().get(index[i]).deviceId.get(j);
					init.getPopList().get(index[i]).deviceId.set(j,
							init.getPopList().get(index[i + 1]).deviceId.get(j));
					init.getPopList().get(index[i + 1]).deviceId.set(j, temp);
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
		int index1 = 0; // 生成一个随机数用于染色体上的某一个点发生变异
		int index2 = 0; // 生成一个随机数，在服务集合里面随机找一个不同的设备，使得该设备发生变化
		for (int i = 0; i < init.getPopList().size(); i++) {// 遍历染色体
			pro = random.nextDouble();
			if (pro < 0.05) {
				index1 = random.nextInt(init.getNumService());// 产生的随机数的范围是：服务数量-1
				deviceId1 = init.getPopList().get(i).deviceId.get(index1);// 记录下来此时的设备Id，从而获得该设备ID的对应服务。
				serviceId1 = init.serviceReq.get(index1);// 记录下来此时的服务Id，从该服务类的集合中找设备
				// 从该服务Id中随机产生一个设备号
				// 若随机产生的设备号和该设备号不同时，替换他，否则，重新产生随机数。
				index2 = random.nextInt(init.getService().get(serviceId1).getDevice().size());
				deviceId2 = init.getService().get(serviceId1).getDevice().get(index2).getId();
				while (deviceId1 == deviceId2) {
					index2 = random.nextInt(init.getService().get(serviceId1).getDevice().size());
					deviceId2 = init.getService().get(serviceId1).getDevice().get(index2).getId();
				}
				init.getPopList().get(i).deviceId.set(index1, deviceId2);
			}
		}
//			System.out.println("测试开始");
//			for(int i = 0; i < init.getPopList().size(); i++ ) {
//				for (int j = 0; j < init.getNumService(); j++) {
//					System.out.print(init.getPopList().get(i).deviceId.get(j) + " ");
//				}
//				System.out.println();
//			}
	}

	/*
	 * 运行算法逻辑,图一二三
	 */
	/*public ArrayList<Double> runNSGA2() {
		// 初始化种群 //贪婪、随机、Nsgaii迁移和不迁移
		// calinsideFitness();
		// init.initPopulation();
		ArrayList<Double> U_max3 = new ArrayList<>();
		// 计算当前种群的适应度
		calOutsideFitness();
		// 将种群保存起来供下次迭代
		copyPopListToPrePopList();
		int gen =0;
		while (gen < maxGen) {
			// 进入下一代
			gen += 1;
			// 交叉变异操作
			crossoverOperator();
			mutationOperator();
			// 计算目前适应度
			calOutsideFitness();
			// 找出支配高，拥挤距离大的个体形成新的种群
			// System.out.println("第"+gen+"代");
			makeNewPopulation();
			//前三幅图需要用到的循环
			if (gen >= 10) {
				System.out.println("第" + gen + "贪婪迁移外部最优");
				double max_u = 0.0;
				max_u = init.chooseOutside();
				U_max3.add(max_u);
			}
			//前三幅图需要用到的循环
			// 将选取出来的种群保存到上一代中
			copyPopListToPrePopList();
		}
		return U_max3;
	}*/
	/*
	 * 图四和图五、图八和图九需要用到的算法运行逻辑
	 */
	public Map runNSGA2(){
		Map greedyResult = new HashMap();

		double max_u= 0.0;
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
			//将选取出来的种群保存到上一代中
			copyPopListToPrePopList();
		}
		greedyResult = init.chooseOutside();
		//max_u = init.chooseOutside();
		return greedyResult;
	}
//	public static void main(String[] args) { 
//		GreedyMig G1 = new GreedyMig(init,50,0.65,0.05,network);
//		
//	}
	 
}