package nsgaii;

import java.util.*;

import nsgaii.Chromsome;
import networkInit.CreateNetwork;
import networkInit.DeviceInit;
import networkInit.Round;
import networkInit.ServiceInit;
import variable.Variable;


public class InitPopulation {
	Random random = new Random();
	CreateNetwork network;
	private int popSize;
	private int numService;
	private int lowBand;
	private int highBand;
	//Round iteRound;
	Vector<Integer> serviceReq;
	private ArrayList<ServiceInit> service;
	private ArrayList<DeviceInit> device;
	private ArrayList<Chromsome> popList;// 当代种群
	private ArrayList<Chromsome> prePopList;// 上一代种群
	HashMap<Chromsome,ArrayList<Chromsome>> outsideDominatingMap;
	public static ArrayList<ArrayList> communicationRange = new ArrayList<>();
	public static ArrayList serviceFlag = new ArrayList();


	public int getNumService() {
		return numService;
	}

	public void setNumService(int numService) {
		this.numService = numService;
	}

	public Vector<Integer> getServiceReq() {
		return serviceReq;
	}

	public void setServiceReq(Vector<Integer> serviceReq) {
		this.serviceReq = serviceReq;
	}

	public InitPopulation(int popSize,int numService,CreateNetwork network,Round iteRound,int lowBand,int highBand) {
		this.popSize = popSize;
		this.numService = numService;
		this.network = network;
		serviceReq = new Vector<Integer>();
		//this.iteRound = iteRound;
		this.lowBand = lowBand;
		this.highBand = highBand;
	}

	public ArrayList<Chromsome> getPopList() {
		return popList;
	}

	public void setPopList(ArrayList<Chromsome> popList) {
		this.popList = popList;
	}

	public ArrayList<Chromsome> getPrePopList() {
		return prePopList;
	}

	public void setPrePopList(ArrayList<Chromsome> prePopList) {
		this.prePopList = prePopList;
	}

	public ArrayList<DeviceInit> getDevice() {
		return device;
	}

	public void setDevice(ArrayList<DeviceInit> device) {
		this.device = device;
	}

	public ArrayList<ServiceInit> getService() {
		return service;
	}

	public void setService(ArrayList<ServiceInit> service) {
		this.service = service;
	}

	public void initPopulation() {
		this.service = new ArrayList<ServiceInit>();
		this.device = new ArrayList<DeviceInit>();
		//int[] serviceRq = new int[numSerivce];
		//network = new CreateNetwork(numAllDevice, numAllService, serviceSkewness, deviceSkewness);
		System.out.println("numService = "+numService);
		network.initService(this.service);//标记几个service->结构框架
		network.initDevice(this.device);//标记几个device->结构框架
		network.allocate(this.service, this.device);//循环所有device使其每个都有service任务
		network.generateSkewnessService(service);//分成计算密集型和数据密集型进行初始化 flag inv wkd cr dt
		network.generateSkewnessDevice(device);//分成普通设备和精英设备进行初始化 x y r flag rsd stg f

		this.popList = new ArrayList<Chromsome>();
		for(int i = 0;i < this.popSize;i++) {
			Chromsome chromsome = new Chromsome(10);
			popList.add(chromsome);
		}
		this.prePopList = new ArrayList<Chromsome>();
		for(int i = 0;i < this.popSize;i++) {
			Chromsome chromsome = new Chromsome(10);
			prePopList.add(chromsome);
		}
		
		float locRet = 0;
		for(int i = 0;i < service.size();i++) {
			for(int j = service.get(i).getDevice().size() - 1;j >= 0;j--) {
			    //service-network distance
				float d = (float) Math.sqrt((network.getIteRound().getX() - service.get(i).getDevice().get(j).getX()) * (network.getIteRound().getX() - service.get(i).getDevice().get(j).getX())
						+ (network.getIteRound().getY() - service.get(i).getDevice().get(j).getY()) * (network.getIteRound().getY() - service.get(i).getDevice().get(j).getY()));
				if (d >= network.getIteRound().getRadius() + Variable.r)
					locRet = 0;
				else if (d <= Math.abs(network.getIteRound().getRadius() - Variable.r)) {
					float r = network.getIteRound().getRadius() < Variable.r ? network.getIteRound().getRadius() : Variable.r;
					locRet = (float) (Math.PI * r * r);
				}
				else{
					float ang1 = (float) Math.acos((network.getIteRound().getRadius() * network.getIteRound().getRadius() + d * d - Variable.r * Variable.r)
							/ 2. / network.getIteRound().getRadius() / d);
					float ang2 = (float) Math.acos(
							(Variable.r * Variable.r + d * d - network.getIteRound().getRadius() * network.getIteRound().getRadius()) / 2. / Variable.r / d);
					float ret = (float) (ang1 * network.getIteRound().getRadius() * network.getIteRound().getRadius() + ang2 * Variable.r * Variable.r
							- d * network.getIteRound().getRadius() * Math.sin(ang1));
					locRet += (float) (ret / (Math.PI * Variable.r * Variable.r));
				}
				if(locRet == 0) {
					service.get(i).getDevice().remove(j);
			    }
			}
		}
		System.out.println("---------service getDevice------------");
		communicationRange.clear();
		serviceFlag.clear();
		for(int i = 0; i < service.size(); i++) {
			//get service flag
			serviceFlag.add(service.get(i).getFlag());
			//get communicationRange
			ArrayList currentDeviceList = new ArrayList();
			communicationRange.add(currentDeviceList);
			for(int j = 0;j < service.get(i).getDevice().size();j++) {
				System.out.print(service.get(i).getDevice().get(j).getId()+" ");
				communicationRange.get(i).add(service.get(i).getDevice().get(j).getId());
			}
			System.out.println();
		}
		int maxCCTRS = 0;
		ArrayList maxCCTR = new ArrayList();
		for(int i=0;i<communicationRange.size();i++){
			if(communicationRange.get(i).size()>maxCCTRS) {
				maxCCTRS = communicationRange.get(i).size();
				maxCCTR = communicationRange.get(i);
			}
		}
		int size = communicationRange.size();
		communicationRange.clear();
		for(int i=0;i<size;i++){
			communicationRange.add(maxCCTR);
		}

		System.out.println();
		for(int i=0;i<service.size();i++){
			for(int j=0;j<service.get(i).getDevice().size();j++){
				service.get(i).getDevice().get(j).setId((Integer) maxCCTR.get(j));
			}

		}






		//产生不同的随机数,作为用户请求的服务的id
		//Vector<Integer> serviceReq = new Vector<Integer>();
		int count = 0;

		// can't get out of this loop!!!!!!!!!
		// TUT
		while(count < numService) {
			int num = random.nextInt(numService);
			//ServiceInit serviceInit = new ServiceInit(serviceReq.get(i),new ArrayList<DeviceInit>());
			if(!serviceReq.contains(num) && service.get(num).getDevice().size()!=0) {
				serviceReq.add(num);
				count++;
			}
			else {
				 continue;
			}
		}

		System.out.println("Initial population:");
		for(int i = 0;i < popList.size();i++) {
			for(int j = 0;j < serviceReq.size();j++) {
				int number = service.get(serviceReq.get(j)).getDevice().size();
//				for (int k = 0;k<number;k++) {
//					System.out.println(service.get(serviceReq.get(j)).getDevice().get(random.nextInt(number)).getId()+" ");
//				}
				popList.get(i).deviceId.add(service.get(serviceReq.get(j)).getDevice().get(random.nextInt(number)).getId()) ;
				popList.get(i).mgDeviceId.add(0);
				System.out.print(popList.get(i).deviceId.get(j)+" ");
			}
			System.out.println();
		}
	}
	
	public void copyToPre() {
		for(int i = 0;i < popSize;i++) {
			prePopList.set(i, (Chromsome)popList.get(i).clone());
		}
	}
	
	public void crossover() {
		Random random = new Random();
		int point;// 配对过程的中间变量
		int temp;// 交换过程中的辅助变量
		
		// 对种群中的个体进行随机配对
		int[] index = new int[popSize];
		for (int i = 0; i < popSize; i++) {
			index[i] = i;
		}
		for (int i = 0; i < popSize; i++) {
			point = random.nextInt(popSize - i);
			temp = index[i];
			index[i] = index[point + i];
			index[point + i] = temp;
		}
		for (int i = 0; i < popSize - 1; i += 2) {
			double pro = random.nextDouble();
			if (pro < 0.65){
				point = random.nextInt(numService);// 交叉点
				// 交换index[i]和index[i+1]在交叉点后面的个体
				for (int j = point; j < numService; j++) {
					temp = popList.get(index[i]).deviceId.get(j);
					popList.get(index[i]).deviceId.set(j, popList.get(index[i + 1]).deviceId.get(j));
					popList.get(index[i + 1]).deviceId.set(j,temp);
				}
			}
		}
		
	}
	
	public void mutation() {
		Random random = new Random();
		Random random2 = new Random();
		double pro = 0.0;// 生成一个随机数用于根据变异概率pm判断当前个体是否参与变异
		int at = 0;
		for (int i = 0; i < popSize; i++) {
			at = random.nextInt(numService);
			pro = random.nextDouble();
			if (pro < 0.05) {
				int a = service.get(serviceReq.get(at)).getDevice().size();
				popList.get(i).deviceId.set(at,
						service.get(serviceReq.get(at)).getDevice().get(random2.nextInt(a)).getId());
			}
		}

	}
	
	public int isDominate(Chromsome p, Chromsome q) {
		/*p = new Chromsome(3);
		q = new Chromsome(3);*/
		int tag = 0;// tag=0表示p和q之间等价，没有支配关系
		if(p.fitness[0] < q.fitness[0] && p.fitness[1] < q.fitness[1] && p.fitness[2] > q.fitness[2] && p.fitness[3] > q.fitness[3]) {
			tag = 1;//p支配q
		}
		else if(p.fitness[0] > q.fitness[0] && p.fitness[1] > q.fitness[1] && p.fitness[2] < q.fitness[2] && p.fitness[3] < q.fitness[3]){
			tag = 2;//q支配p
		}
		return tag;
	}
	
	public ArrayList<Chromsome> getSortList(ArrayList<Chromsome> list,int numObjective) {
		Chromsome temp = new Chromsome();
		for (int i = 0; i < list.size(); i++) {
			for (int j = i + 1; j < list.size(); j++) {
				if (list.get(i).fitness[numObjective] > list.get(j).fitness[numObjective]) {
					temp = list.get(i);
					list.set(i, (Chromsome)list.get(j).clone());
					list.set(j, (Chromsome)temp.clone());
				}
			}
		}
		return list;
	}
	
	public ArrayList<ArrayList<Chromsome>> fastNondominate() {

		outsideDominatingMap = new HashMap<Chromsome,ArrayList<Chromsome>>();
		// 第一个pareto等级集合
		ArrayList<Chromsome> firstParetoRankSet = new ArrayList<Chromsome>();
		// 所有pareto等级集合
		ArrayList<ArrayList<Chromsome>> paretoRankSetList = new ArrayList<ArrayList<Chromsome>>();// pareto分级的解集
		// 两代种群集合
		ArrayList<Chromsome> unionPopList = new ArrayList<Chromsome>();
		for (int i = 0; i < 2 * popSize; i++) {
			Chromsome chromsome = new Chromsome(4);
			unionPopList.add(chromsome);
		}
/*		for (int i = 0; i < popSize; i++) {
		unionPopList.set(i, prePopList.get(i));
		}
		for (int i = popSize; i < 2*popSize; i++ ){
			unionPopList.set(i, popList.get(i - popSize));
		}
*/
		// 将prePopList中的元素复制到unionPopList中
		for (int i = 0; i < popSize; i++) {
			unionPopList.set(i, (Chromsome)prePopList.get(i).clone());
		}
			
		// 将PopList中的元素复制到unionPopList中
		for (int i = popSize; i < (2 * popSize); i++) {
			unionPopList.set(i,(Chromsome)popList.get(i - popSize).clone());
		}
		// 每个个体和所有个体比较,判断支配
		for (int i = 0; i < unionPopList.size(); i++) {
			ArrayList<Chromsome> tempDominatingList = new ArrayList<Chromsome>();
			for (int j = 0; j < unionPopList.size(); j++) {
				int tag = this.isDominate(unionPopList.get(i), unionPopList.get(j));
				if (tag == 1) {
					//unionPopList.get(i).getDominatingList().add(unionPopList.get(j));// 记录当前个体支配的解
					tempDominatingList.add(unionPopList.get(j));
				} else if (tag == 2) {
					unionPopList.get(i).numDominated += 1;// 支配当前个体的个体数目加1
				}
			}
			outsideDominatingMap.put(unionPopList.get(i), tempDominatingList);
			// 如果当前解的numDominated属性等于0，即属于支配最前沿，则将其加入第一个pareto等级集合，并将其paretoRank值设置为1
			if (unionPopList.get(i).numDominated == 0) {
				unionPopList.get(i).setParetoRank(1);
				firstParetoRankSet.add(unionPopList.get(i));
			}
		}
		paretoRankSetList.add(firstParetoRankSet);
		//System.out.println("pareto front 集合的大小： " + firstParetoRankSet.size());

		int rank = 0;
		while (paretoRankSetList.get(rank).size() > 0
				&& paretoRankSetList.get(rank) != null) {
			// 用于储存下一前沿的pareto集合
			ArrayList<Chromsome> paretoRankSet = new ArrayList<Chromsome>();
			// 依次处理当前pareto等级里的所有个体
			for (int j = 0; j < paretoRankSetList.get(rank).size(); j++) {
				//Chromsome currentChromsome = paretoRankSetList.get(rank).get(j);
				ArrayList<Chromsome> currentDList = outsideDominatingMap.get(paretoRankSetList.get(rank).get(j));
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
		return paretoRankSetList;
	}

	public void outsideCrowdingDistance(ArrayList<ArrayList<Chromsome>> paretoRankSetList) {
		for (int i = 0; i < paretoRankSetList.size(); i++) {
			ArrayList<Chromsome> paretoRankSet = new ArrayList<Chromsome>();
			paretoRankSet = paretoRankSetList.get(i);

			if (paretoRankSet.size() > 0 && paretoRankSet != null) {
				//按照目标依次计算个体的crowdingDistance值
				for (int j = 0; j < 4; j++) {
					//对paretoRankSet中的个体按照fitness从小到大进行排序
					paretoRankSet = this.getSortList(paretoRankSet, j);
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
				Chromsome temp = new Chromsome(4);
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
	
	public void makeOutsideNewPopulation() {
		// 对两代种群快速非支配排序
		ArrayList<ArrayList<Chromsome>> paretoRankSetList = this.fastNondominate();
		// 对每个pareto前沿计算拥挤距离并排序
		this.outsideCrowdingDistance(paretoRankSetList);
		// 记录加入到popList中的个体总数
		int numChromsome = 0;
		// popList的下标
		int popListIndex = 0;
		//取出前面n个个体作为下一代种群
		for (int i = 0; i < paretoRankSetList.size(); i++) {
			ArrayList<Chromsome> paretoRankSet = paretoRankSetList.get(i);
			numChromsome += paretoRankSet.size();

			if (numChromsome < popSize && popListIndex < popSize) {
				for (int j = 0; j < paretoRankSet.size(); j++) {
					popList.set(popListIndex++, paretoRankSet.get(j));
					// popList.get(popListIndex++).setBinChrom(paretoRankSet.get(j).getBinChrom());
				}
			} else {
				int overNum = numChromsome - popSize;
				for (int j = 0; j < paretoRankSet.size() - overNum
						&& popListIndex < popSize; j++) {
					popList.set(popListIndex++, paretoRankSet.get(j));
					// popList.get(popListIndex++).setBinChrom(paretoRankSet.get(j).getBinChrom());
				}
			}
		}
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
/*		//System.out.println("外部最优");
		for(int i = 0;i < bestSolution.size();i++) {
			for(int j = 0;j < bestSolution.get(i).deviceId.size();j++) {
				System.out.print(bestSolution.get(i).deviceId.get(j)+" ");
			}
			System.out.println();
//		System.out.println("对应策略：");
//		System.out.println();
//		for(int i = 0;i < bestSolution.size();i++) {
//			for(int j = 0;j < bestSolution.get(i).mgDeviceId.size();j++) {
//				System.out.print(bestSolution.get(i).mgDeviceId.get(j)+" ");
//			}
//			System.out.println();
		}*/
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
		//画图四和图五时，需要分别输出服务组合能耗conErg和服务组合时间conTime
		double conErg = 0.0;
		double conTime = 0.0;
		//服务迁移时间和用户数量之间的关系，这里定义一个内层的消耗时间变量
		double insideconTime = 0.0;
		double insideconErg = 0.0;
		//服务组合消耗的总能量应是迁移能耗+服务组合能耗，即内层的fitness值和外层的fitness值之和
		//System.out.println("输出内层的适应度函数值： " + bestSolution.get(0).insidefitness[1]);
		//System.out.println("输出外层的适应度函数值： " + bestSolution.get(0).fitness[0]);
		for(int i = 0; i < bestSolution.size(); i++) {
			conErg =  bestSolution.get(i).fitness[0];
			conTime = bestSolution.get(i).fitness[1];
			insideconTime = bestSolution.get(i).insidefitness[0];
			insideconErg = bestSolution.get(i).insidefitness[1];
		}
		//图四、图八
		//System.out.println("conErgde的值是" + conErg);
		//return conErg;
		//图五、图九
		//return conTime;
		//服务迁移时间和用户数量之间的关系
		//return insideconTime;
		//服务迁移时间和用户数量之间的关系


		Map result = new HashMap();

		result.put("bestSolution",bestSolution);
		result.put("conErg",conErg);
		result.put("conTime",conTime);
		return result;
		//return insideconErg;
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
	
	public float calDij(DeviceInit device1,DeviceInit device2) {
		//距离
		return (float) Math.sqrt((device1.getX() - device2.getX()) * (device1.getX() - device2.getX()) + (device1.getY() - device2.getY()) * (device1.getY() - device2.getY()));
	}
	
	public float calNoise(DeviceInit device1,DeviceInit device2) {
		//噪声功率
		if(device1.getFlag() == 1 && device2.getFlag() == 1) {
			return (float) (Math.pow(10, ((-174 + 10 * (Math.log(lowBand) / Math.log(10))) / 10)) / 1000);
		}
		else
		return (float) (Math.pow(10, ((-174 + 10 * (Math.log(highBand) / Math.log(10))) /10 )) / 1000);
	}
	
	public float calRij(DeviceInit device1,DeviceInit device2) {
		//rij
		float rij = 0.0f;
		if(device1.getFlag() == 1 && device2.getFlag() == 1) {
			rij = (float) (lowBand * Math.log(1 + 0.2 * Math.pow(calDij(device1, device2), -4) / calNoise(device1, device2)) / Math.log(2));
		}
		else {
			rij = (float) (highBand * Math.log(1 + 0.2 * Math.pow(calDij(device1, device2), -4) / calNoise(device1, device2)) / Math.log(2));
		}
		return rij;
	}
	
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
	
	public void calNoMigrationFitness(){
		for(int i=0;i<popList.size();i++) {
			for(int j=0;j<popList.get(i).deviceId.size();j++) {
				popList.get(i).mgDeviceId.set(j, popList.get(i).deviceId.get(j));
			}
		}
		for(int i=0;i<popList.size();i++) {
			popList.get(i).insidefitness[0] = 0.0;
			popList.get(i).insidefitness[1] = 0.0;
		}
		for(int i = 0;i < popList.size();i++) {
			double e = 0.0;
			double e_inv = 0.0;
			double e_cmp = 0.0;
			double e_trs = 0.0;
			double t = 0.0;
			double t_cmp = 0.0;
			double t_trs = 0.0;
			double spt = 0.0;
			double lbn = 0.0;
			double lbd = 0.0;
			for(int j = 0;j < popList.get(i).deviceId.size();j++) {
				e_inv += service.get(serviceReq.get(j)).getInv();
				e_cmp += Math.pow(10, -26) * device.get(popList.get(i).deviceId.get(j)).getF() * device.get(popList.get(i).deviceId.get(j)).getF()
						* service.get(serviceReq.get(j)).getCr();
				t_cmp += service.get(serviceReq.get(j)).getCr()
						/ device.get(popList.get(i).deviceId.get(j)).getF();
				spt += evaluateLoc(device.get(popList.get(i).deviceId.get(j)), network.getIteRound());
				
				if(j != (popList.get(i).deviceId.size() - 1)) {
					e_trs += (0.2 * service.get(serviceReq.get(j)).getDt() 
							/ calRij(device.get(popList.get(i).deviceId.get(j)), device.get(popList.get(i).deviceId.get(j + 1))));
					t_trs += service.get(serviceReq.get(j)).getDt() 
							/ calRij(device.get(popList.get(i).deviceId.get(j)), device.get(popList.get(i).deviceId.get(j + 1)));
				}
			}
			e += e_inv + e_cmp + e_trs;
			t = t_cmp + t_trs;
			spt = spt / popList.get(i).deviceId.size();
			for(int m = 0;m < popList.get(i).deviceId.size();m++) {
				lbd += (device.get(popList.get(i).deviceId.get(m)).getRsd() - e) / device.get(popList.get(i).deviceId.get(m)).getRsd();
			}
			lbn = lbd / popList.get(i).deviceId.size();
			popList.get(i).fitness[0] = e;
			popList.get(i).fitness[1] = t;
			popList.get(i).fitness[2] = spt;
			popList.get(i).fitness[3] = lbn;
		}

	}
	/*public static void main(String[] args) {
		InitPopulation init = new InitPopulation(20);
		//ArrayList<ServiceInit> service = new ArrayList<ServiceInit>();
		//ArrayList<DeviceInit> device = new ArrayList<DeviceInit>();
		//ArrayList<Chromsome> popList = new ArrayList<Chromsome>();// 当代种群
		init.initPopulation();
	//	System.out.print(popList.get(i).deviceId[j]+"");
	}*/
}
