package nsgaii;

import java.io.IOException;
import java.util.*;

import networkInit.Round;
import networkInit.CreateNetwork;
import networkInit.DeviceInit;
import variable.Variable;
import nsgaii.InitPopulation;


public class Inside  {
	   
	
	private int insidePopSize;
	private int maxGen;
	private int lowBand;
	private int highBand;
	ArrayList<ArrayList<Chromsome>> migrationPopList;
	ArrayList<ArrayList<Chromsome>> migrationPrePopList;
	ArrayList<Chromsome> insidePopList;
	//ArrayList<ArrayList<ArrayList<Chromsome>>> paretoRankSetListAll = new ArrayList<ArrayList<ArrayList<Chromsome>>>();
	InitPopulation init;
	
	private int numService;
	CreateNetwork network;
	HashMap<Chromsome,ArrayList<Chromsome>> dominatingMap;

	private int[][] chain;
	private int[] migration;
	
	public Inside(int insidePopSize, InitPopulation init, int numService, CreateNetwork network,int maxGen,int lowBand,int highBand) {
		this.insidePopSize = insidePopSize;
		this.init = init;
		this.numService = numService;
		this.network = network;
		this.maxGen = maxGen;
		this.lowBand = lowBand;
		this.highBand = highBand;
	}
	
	public void migrationMatrix(int[][] chains) {
		Random random = new Random();
		for(int j = 0;j < numService;j++) {
			int a = random.nextInt(numService);
			for(int i = 0;i < numService;i++) {
				if(i == a) {
					chains[i][j] = 1;
				}
				else {
					chains[i][j] = 0;
				}
				
			}
		}
		/*for(int i = 0;i < b;i++) {
			for(int j = 0;j < b;j++) {
				System.out.print(chains[i][j]+" ");
			}
			System.out.println();
		}*/
	}

	public int[] calMigration(ArrayList<Integer> x,int[][] y) {
		int c[] = new int[numService];
		for(int i = 0;i < numService;i++) {
			for(int j = 0;j < numService;j++) {
				c[i] += x.get(j)*y[j][i]; 
			}
		}
		return c;
	}
	
	public void initInsidePopulation() {
		chain = new int[numService][numService];
		migration = new int[numService];
		migrationPopList = new ArrayList<ArrayList<Chromsome>>();

		for(int i = 0;i < numService;i++) {
			migration[i] = 0;
		}

		for(int i = 0;i < init.getPopList().size();i++) {
			insidePopList = new ArrayList<Chromsome>();
			
			for(int j = 0;j < insidePopSize;j++) {
				Chromsome chromsome = new Chromsome(10);
				insidePopList.add(chromsome);
				migrationMatrix(chain);
				migration = calMigration(init.getPopList().get(i).deviceId,chain).clone();
				for(int l = 0;l < migration.length;l++) {
					insidePopList.get(j).deviceId.add(migration[l]);
				}
			}
			migrationPopList.add(insidePopList);
		}
	}
	
	public float calDij(DeviceInit device1,DeviceInit device2) {
		////距离
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
	
	public void evaluateInside() {
		int count = 0;
		float c_cst = 0.0f;
		float e_trs = 0.0f;
		for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < insidePopSize;j++) {
				for(int m = 0;m < numService;m++) {
					//if(init.getPopList().get(i).deviceId.get(m) != migrationPopList.get(i).get(j).deviceId.get(m)) {
					for(int n = 0;n < numService;n++) {
						if(migrationPopList.get(i).get(j).deviceId.get(m).equals(migrationPopList.get(i).get(j).deviceId.get(n))) {
							c_cst += init.getService().get(init.getServiceReq().get(n)).getWkd();
							count++;
						}
						if(migrationPopList.get(i).get(j).deviceId.get(m).equals(migrationPopList.get(i).get(j).deviceId.get(n)) &&
							!init.getPopList().get(i).deviceId.get(n).equals(migrationPopList.get(i).get(j).deviceId.get(n))) {
							//e_inv += init.getService().get(init.getServiceReq().get(n)).getInv();
							//1g = 1048576 kb
							e_trs += 0.2 * init.getService().get(init.getServiceReq().get(n)).getWkd() * 1048576 / calRij(init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(n)), init.getDevice().get(init.getPopList().get(i).deviceId.get(n)));
							//System.out.println("Rij: "+calRij(init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(n)), init.getDevice().get(init.getPopList().get(i).deviceId.get(n))));
						}
					}
					while(c_cst > init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m)).getStg() &&
							count > 5 &&
							e_trs > init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m)).getRsd()) {
						migrationMatrix(chain);
						migration = calMigration(init.getPopList().get(i).deviceId,chain).clone();
						for(int l = 0;l < migration.length;l++) {
							migrationPopList.get(i).get(j).deviceId.set(l,migration[l]);
						}
						count = 0;
						c_cst = 0.0f;
						e_trs = 0.0f;
						for(int n = 0;n < numService;n++) {
							if(migrationPopList.get(i).get(j).deviceId.get(m).equals(migrationPopList.get(i).get(j).deviceId.get(n))) {
								c_cst += init.getService().get(init.getServiceReq().get(n)).getWkd();
								count++;
							}
							if(migrationPopList.get(i).get(j).deviceId.get(m).equals(migrationPopList.get(i).get(j).deviceId.get(n))  &&
								!init.getPopList().get(i).deviceId.get(n).equals(migrationPopList.get(i).get(j).deviceId.get(n)) ) {
								e_trs += 0.2 * init.getService().get(init.getServiceReq().get(n)).getWkd() / calRij(init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(n)), init.getDevice().get(init.getPopList().get(i).deviceId.get(n)));

							}
						}
					}
						
					count = 0;
					c_cst = 0.0f;
					e_trs = 0.0f;
				}
					
			}
		}
	//}
	/*	System.out.println("初始");
		for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < 10;j++) {
				for(int l = 0;l < numService;l++) {
					System.out.print(migrationPopList.get(i).get(j).deviceId.get(l)+" ");
				}
				System.out.println();
			}
			System.out.println();
		}
		*/
}
	
	public void calFitness() {
		//evaluateInside();
		double t = 0.0;
		double e = 0.0;
		double cst = 0.0;
		double cbf = 0.0;
		//int count = 0;
		for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < insidePopSize;j++) {
				ArrayList<Integer> temp = new ArrayList<Integer>();
				ArrayList<Double> temp_cst = new ArrayList<Double>();
				//double[] fitness = new double[3];
				for(int m = 0;m < numService;m++) {
					for(int n = 0;n < numService;n++) {
						if(migrationPopList.get(i).get(j).deviceId.get(m).equals(migrationPopList.get(i).get(j).deviceId.get(n)) 
								&& !temp.contains(migrationPopList.get(i).get(j).deviceId.get(m))) {
							cst += init.getService().get(init.getServiceReq().get(n)).getWkd();	
							//count++;
						}
					}
					if(!temp.contains(migrationPopList.get(i).get(j).deviceId.get(m))) {
						temp_cst.add(cst);
						temp.add(migrationPopList.get(i).get(j).deviceId.get(m));
					}
					cst = 0.0;
					if(!init.getPopList().get(i).deviceId.get(m).equals(migrationPopList.get(i).get(j).deviceId.get(m))) {
						t += init.getService().get(init.getServiceReq().get(m)).getWkd() * 1048576 / calRij(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)), init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m)))
								+ init.getService().get(init.getServiceReq().get(m)).getCr() / init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m)).getF();
						/*System.out.println("距离"+calDij(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)), init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m))));
						System.out.println("噪声"+calNoise(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)), init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m))));
						System.out.println(calRij(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)), init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m))));*/
						e += 0.2 * init.getService().get(init.getServiceReq().get(m)).getWkd() * 1048576/ calRij(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)), init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m)))
								+ 0.002 * init.getService().get(init.getServiceReq().get(m)).getCr() / init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m)).getF();
					}
					else {
						if(m < (numService - 1)) {
							t += init.getService().get(init.getServiceReq().get(m)).getCr() / init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m)).getF()
									+ init.getService().get(init.getServiceReq().get(m)).getDt() / calRij(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)), init.getDevice().get(init.getPopList().get(i).deviceId.get(m+1)));
							e += 0.2 * init.getService().get(init.getServiceReq().get(m)).getDt() / calRij(init.getDevice().get(init.getPopList().get(i).deviceId.get(m)), init.getDevice().get(init.getPopList().get(i).deviceId.get(m+1)))
									+ 0.5 * init.getService().get(init.getServiceReq().get(m)).getCr() / init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m)).getF();
						}
						else if(m == (numService - 1)){
							t += init.getService().get(init.getServiceReq().get(m)).getCr() / init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m)).getF();
							e += 0.5 * init.getService().get(init.getServiceReq().get(m)).getCr() / init.getDevice().get(migrationPopList.get(i).get(j).deviceId.get(m)).getF();
						}
						
					}
				}
				for(int x = 0;x < temp.size();x++) {
					cbf += (init.getDevice().get(temp.get(x)).getRsd() - temp_cst.get(x)) / init.getDevice().get(temp.get(x)).getRsd();
				
				}
				cbf = cbf / numService;
				migrationPopList.get(i).get(j).fitness[0] = t;
				migrationPopList.get(i).get(j).fitness[1] = e;
				migrationPopList.get(i).get(j).fitness[2] = cbf;

				t = 0.0;
				e = 0.0;
				cbf = 0.0;
			}
		}
		/*System.out.println("内层的fitness");
		for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < 10;j++) {
				for(int m = 0;m < 3;m++) {
					System.out.print(migrationPopList.get(i).get(j).fitness[m]+" ");
				}
				System.out.println();
			}
			System.out.println();
		}*/
	}

	public void copyPopListToPrePopList() {
		this.migrationPrePopList = new ArrayList<ArrayList<Chromsome>>();
		for(int i = 0;i < init.getPopList().size();i++) {
			insidePopList = new ArrayList<Chromsome>();
			for(int j = 0;j < insidePopSize;j++) {
				Chromsome chromsome = new Chromsome(3);
				insidePopList.add(chromsome);
			}
			migrationPrePopList.add(insidePopList);
		}
		for (int i = 0; i < init.getPopList().size(); i++) {
			for(int j = 0;j < insidePopSize;j++) {
				//migrationPopList.get(i).get(j).numDominated = 1;
				//System.out.println(System.identityHashCode(migrationPrePopList.get(i).get(j).dominatingList));
				//System.out.println("dizhi1="+System.identityHashCode(migrationPopList.get(i).get(j).dominatingList));
				//migrationPrePopList.get(i).set(j, migrationPopList.get(i).get(j));
				migrationPrePopList.get(i).set(j, (Chromsome)migrationPopList.get(i).get(j).clone());
				//System.out.println("dizhi2="+System.identityHashCode(migrationPrePopList.get(i).get(j).dominatingList));
			}
		}
		/*System.out.println("上一代");
		for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < 10;j++) {
				for(int m = 0;m < b;m++) {
					System.out.print(migrationPrePopList.get(i).get(j).deviceId.get(m)+" ");
					//System.out.print("num"+migrationPrePopList.get(i).get(j).paretoRank+" ");
				}
				System.out.println();
			}
			System.out.println();
		}*/
	}
	
	public void crossoverOperator() {
		Random random = new Random();
		int point;// 配对过程的中间变量
		int temp;// 交换过程中的辅助变量

		// 对种群中的个体进行随机配对
		int[] index = new int[insidePopSize];
		for (int i = 0; i < insidePopSize; i++) {
			index[i] = i;
		}
		for (int i = 0; i < insidePopSize; i++) {
			point = random.nextInt(insidePopSize - i);
			temp = index[i];
			index[i] = index[point + i];
			index[point + i] = temp;
		}
		for(int i = 0;i < init.getPopList().size();i++) {
			for (int j = 0; j < insidePopSize - 1; j += 2) {
				double pro = random.nextDouble();
				if (pro < 0.65){
					point = random.nextInt(numService);// 交叉点
					// 交换index[i]和index[i+1]在交叉点后面的个体
					for (int m = point; m < numService; m++) {
						temp = migrationPopList.get(i).get(index[j]).deviceId.get(m);
						migrationPopList.get(i).get(index[j]).deviceId.set(m,migrationPopList.get(i).get(index[j + 1]).deviceId.get(m));
						migrationPopList.get(i).get(index[j + 1]).deviceId.set(m,temp);
					}
				}
			}
		}
		/*System.out.println("交叉后");
		for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < 10;j++) {
				for(int m = 0;m < numService;m++) {
					System.out.print(migrationPopList.get(i).get(j).deviceId.get(m)+" ");
				}
				System.out.println();
			}
			System.out.println();
		}*/
		/*System.out.println("交叉后pre");
		for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < 10;j++) {
				for(int m = 0;m < b;m++) {
					System.out.print(migrationPrePopList.get(i).get(j).deviceId.get(m)+" ");
				}
				System.out.println();
			}
			System.out.println();
		}*/
	}

	// 按概率pm进行单点变异，即随机确定一个变异点然后取反
	public void mutationOperator() {
		Random random = new Random();
		//Random random2 = new Random();
		double pro = 0.0;// 生成一个随机数用于根据变异概率pm判断当前个体是否参与变异
		for(int i = 0;i < init.getPopList().size();i++) {
			for (int j = 0; j < insidePopSize; j++) {
				for (int m = 0; m < numService; m++) {
					pro = random.nextDouble();
					if (pro < 0.05) {
						migrationPopList.get(i).get(j).deviceId.set(m,migrationPopList.get(i).get(j).deviceId.get(random.nextInt(numService)));
					}
				}
			}
		}
	}

	public void repair() {
		boolean overfilled = false;
		Random random = new Random();
		for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < insidePopSize;j++) {
				for(int m = 0;m < numService;m++) {
					int count = 0;
					int id = 0;
					ArrayList<Integer> temp = new ArrayList<>();
					ArrayList<Integer> temp2 = new ArrayList<>();
					for(int n = 0;n < numService;n++) {
						if(migrationPopList.get(i).get(j).deviceId.get(n).equals(migrationPopList.get(i).get(j).deviceId.get(m))) {
							temp.add(n);
							count++;
						}
					}
					if(count > 5) {
						overfilled = true;
						id = migrationPopList.get(i).get(j).deviceId.get(m);
					}
					if(overfilled) {
						while(overfilled) {
							for(int z = 0;z < temp.size() - 1;z++) {
								/*if(temp.get(z) + 1 == temp.get(z+1)) {
									temp2.add(temp.get(z));
									temp2.add(temp.get(z+1));

								}*/
								int l = z+1;
								int jump = 1;
								boolean in = false;
								while(l < temp.size() && temp.get(z)+jump == temp.get(l) && !temp2.contains(temp.get(z))) {
									temp2.add(temp.get(l));
									jump++;
									l++;
									in = true;
								}
								if(in) {
									temp2.add(temp.get(z));
								}
							}
							if((temp2.size() != 0)) {
								for(int a = 0;a < temp2.size();a++) {
									if(temp.contains(temp2.get(a))) {
										temp.remove(temp2.get(a));
									}
								}
								if(temp.size() >= temp.size()+temp2.size()-5) {
									for(int b = 0;b < temp.size()+temp2.size()-5;b++) {
										migrationPopList.get(i).get(j).deviceId.set(temp.get(b), migrationPopList.get(i).get(j).deviceId.get(random.nextInt(numService)));
									}
								}
								else if((temp2.size() != 0) && ((temp.size()+temp2.size() - 5) > temp.size())) {
									for(int a = 0;a < temp.size();a++) {
										migrationPopList.get(i).get(j).deviceId.set(temp.get(a), migrationPopList.get(i).get(j).deviceId.get(random.nextInt(numService)));
									}
									for(int b = 0;b < temp2.size()-5;b++) {
										migrationPopList.get(i).get(j).deviceId.set(temp2.get(b), migrationPopList.get(i).get(j).deviceId.get(random.nextInt(numService)));
									}
								}

							}
							else {
								for(int d = 0;d < temp.size() - 5;d++) {
									migrationPopList.get(i).get(j).deviceId.set(temp.get(d), migrationPopList.get(i).get(j).deviceId.get(random.nextInt(numService)));
								}
							}
							int number = 0;
							for(int e = 0;e < numService;e++) {
								if(id == migrationPopList.get(i).get(j).deviceId.get(e)) {
									number++;
								}
							}
							//System.out.println(number);
							//System.out.println(id);
							if(number <= 5) {
								overfilled = false;
								//System.out.println("false");
							}
						}
					}
				}
			}
		}
	}

	public int isDominate(Chromsome p, Chromsome q) {
		/*p = new Chromsome(3);
		q = new Chromsome(3);*/
		int tag = 0;// tag=0表示p和q之间等价，没有支配关系
		if(p.fitness[0] < q.fitness[0] && p.fitness[1] < q.fitness[1] && p.fitness[2] > q.fitness[2]) {
			tag = 1;//p支配q
		}
		else if(p.fitness[0] > q.fitness[0] && p.fitness[1] > q.fitness[1] && p.fitness[2] < q.fitness[2]){
			tag = 2;//q支配p
		}
		return tag;
	}

	 //对list列表中的个体按照fintness值从小到大进行排序(快速排序)
	public ArrayList<Chromsome> getSortList(ArrayList<Chromsome> list,int numObjective) {
		Chromsome temp = new Chromsome(3);
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

	public ArrayList<ArrayList<ArrayList<Chromsome>>> fastNondominateSort() {

		dominatingMap = new HashMap<Chromsome,ArrayList<Chromsome>>();
		// 第一个pareto等级集合
		ArrayList<Chromsome> firstParetoRankSet = new ArrayList<Chromsome>();
		// 所有pareto等级集合
		ArrayList<ArrayList<Chromsome>> paretoRankSetList = new ArrayList<ArrayList<Chromsome>>();// pareto分级的解集
		// 两代种群集合
		ArrayList<ArrayList<Chromsome>> unionPopList = new ArrayList<ArrayList<Chromsome>>();
		ArrayList<ArrayList<ArrayList<Chromsome>>> paretoRankSetListAll = new ArrayList<ArrayList<ArrayList<Chromsome>>>();
		for(int i = 0;i < init.getPopList().size();i++) {
			insidePopList = new ArrayList<Chromsome>();
			for(int j = 0;j < (2 * insidePopSize);j++) {
				Chromsome chromsome = new Chromsome(3);
				insidePopList.add(chromsome);
			}
			unionPopList.add(insidePopList);
		}

		// 将prePopList中的元素复制到unionPopList中
		/*System.out.println("结合前的上一代");
		for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < 10;j++) {
				for(int m = 0;m < b;m++) {
					System.out.print(migrationPrePopList.get(i).get(j).deviceId.get(m)+" ");
					//System.out.print(migrationPrePopList.get(i).get(j).numDominated+" ");
				}
				System.out.println();
			}
			System.out.println();
		}*/
		for(int i = 0;i < init.getPopList().size();i++) {
			for (int j = 0; j < insidePopSize; j++) {
				unionPopList.get(i).set(j,(Chromsome)migrationPrePopList.get(i).get(j).clone());
			}
		}
		/*System.out.println("结合前的当前种群");
		for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < 10;j++) {
				for(int m = 0;m < b;m++) {
					System.out.print(migrationPopList.get(i).get(j).deviceId.get(m)+" ");
				}
				System.out.println();
			}
			System.out.println();
		}*/
		// 将PopList中的元素复制到unionPopList中
		for(int i = 0;i < init.getPopList().size();i++) {
			for (int j = insidePopSize; j < (2 * insidePopSize); j++) {
				//System.out.println("size:"+init.getPopList().size());
				unionPopList.get(i).set(j,(Chromsome)migrationPopList.get(i).get(j - insidePopSize).clone());
				//System.out.println("G");
			}
		}

		/*System.out.println("结合后");
		for(int i = 0;i < init.getPopList().size();i++) {
			for (int j = 0; j < (2 * 10); j++) {
				for (int m = 0; m < b; m++) {
					System.out.print(unionPopList.get(i).get(j).deviceId.get(m) + " ");
					//System.out.print("num"+unionPopList.get(i).get(j).numDominated + " ");
				}
				System.out.println();
				for (int n = 0; n < 3; n++) {
					System.out.print(unionPopList.get(i).get(j).fitness[n] + " ");
				}
				System.out.println();
			}
			System.out.println();
		}*/
		// 每个个体和所有个体比较,判断支配
		for(int i = 0;i < init.getPopList().size();i++) {
			for (int j = 0; j < unionPopList.get(i).size(); j++) {
				ArrayList<Chromsome> tempDominatingList = new ArrayList<Chromsome>();
				//insidePopList = new ArrayList<Chromsome>();
				for (int m = 0; m < unionPopList.get(i).size(); m++) {

					int tag = this.isDominate(unionPopList.get(i).get(j), unionPopList.get(i).get(m));
					if (tag == 1) {
						//unionPopList.get(i).get(j).getDominatingList().add(unionPopList.get(i).get(m));// 记录当前个体支配的解
						tempDominatingList.add(unionPopList.get(i).get(m));

					} else if (tag == 2) {
						unionPopList.get(i).get(j).numDominated += 1;// 支配当前个体的个体数目加1
					}
				}
				dominatingMap.put(unionPopList.get(i).get(j), tempDominatingList);
				// 如果当前解的numDominated属性等于0，即属于支配最前沿，则将其加入第一个pareto等级集合，并将其paretoRank值设置为1
				if (unionPopList.get(i).get(j).numDominated == 0) {
					unionPopList.get(i).get(j).setParetoRank(0);
					//Chromsome chromsome = new Chromsome(3);
					//insidePopList.add(chromsome);
					firstParetoRankSet.add(unionPopList.get(i).get(j));
					//firstParetoRankSet.get(i).add(insidePopList.get(x));
					//x++;
					//paretoRankSetList.add(new ArrayList<Chromsome>());

				}
			}
			paretoRankSetList.add(firstParetoRankSet);
			paretoRankSetListAll.add(paretoRankSetList);
			firstParetoRankSet = new ArrayList<Chromsome>();
			paretoRankSetList = new ArrayList<ArrayList<Chromsome>>();

			//x = 0;
			//System.out.println("GG");
		}

		/*for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < paretoRankSetListAll.get(i).size();j++) {
				for(int m = 0;m < paretoRankSetListAll.get(i).get(j).size();m++) {
					for(int n = 0;n < b;n++) {
						System.out.print(paretoRankSetListAll.get(i).get(j).get(m).deviceId.get(n)+" ");
					}
					System.out.println();
					for(int l = 0;l < 3;l++) {
						System.out.print(paretoRankSetListAll.get(i).get(j).get(m).fitness[l]+" ");
					}
					System.out.println();
				}
			}
			System.out.println();
		}*/


		//System.out.println("pareto front 集合的大小： " + firstParetoRankSet.size());
		//System.out.println(paretoRankSetListAll.get(3).get(0).size());
		int rank = 0;
		for(int i = 0;i < init.getPopList().size();i++) {
			//System.out.println("KKK");
			//for(int m = 0;m < rank;m++) {
				while (paretoRankSetListAll.get(i).get(rank).size() > 0 && paretoRankSetListAll.get(i).get(rank) != null) {
					//for(int m = 0;m < )
					// 用于储存下一前沿的pareto集合
					ArrayList<Chromsome> paretoRankSet = new ArrayList<Chromsome>();
					// 依次处理当前pareto等级里的所有个体
					for (int j = 0; j < paretoRankSetListAll.get(i).get(rank).size(); j++) {
						//Chromsome currentChromsome = paretoRankSetListAll.get(i).get(rank).get(j);
						ArrayList<Chromsome> currentDList = dominatingMap.get(paretoRankSetListAll.get(i).get(rank).get(j));
						if (currentDList.size() > 0 && currentDList != null) {
							for (int k = 0; k < currentDList.size(); k++) {
								Chromsome dominatedChromsome = currentDList.get(k);
								dominatedChromsome.numDominated -= 1;
								// 若当前个体属于下一个前沿，修改个体的paretoRank属性值,并将其加入相应的pareto等级集合
								if (dominatedChromsome.numDominated == 0) {
									dominatedChromsome.setParetoRank(rank + 1);
									paretoRankSet.add(dominatedChromsome);
								}
								//System.out.println("XXXXX");
							}
						}
						//System.out.println("GGGG");
					}
					if (paretoRankSet != null && paretoRankSet.size() > 0) {// 如果paretoRankSet不为空，将其加入paretoRankSetList中；否则，说明所有个体都已经处理完毕，退出循环
						paretoRankSetListAll.get(i).add(paretoRankSet);// 将新生成的paretoRankSet加入paretoRankSetList中供后面计算crowdingDistance时使用
							/*paretoRankSetList.add(paretoRankSet);
						paretoRankSetListAll.add(paretoRankSetList);
						paretoRankSetList = new ArrayList<ArrayList<Chromsome>>();*/
						rank += 1;
					} else {
						//System.out.println("BREAK");
						break;
					}
					//System.out.println("GGG");
				}
				rank = 0;
			//}
				//System.out.println("AAA");
		}

			/*System.out.println("分层后");
			for(int i = 0;i < init.getPopList().size();i++) {
				System.out.println("第" + i + "个：");
				for(int j = 0;j < paretoRankSetListAll.get(i).size();j++) {
					System.out.println("第" + j + "层");
					for(int m = 0;m < paretoRankSetListAll.get(i).get(j).size();m++) {
						for(int n = 0;n < b;n++) {
							System.out.print(paretoRankSetListAll.get(i).get(j).get(m).deviceId.get(n)+" ");
							//System.out.print("num"+paretoRankSetListAll.get(i).get(j).get(m).paretoRank+" ");
						}

						System.out.println();
					}
					System.out.println();
				}
				System.out.println();
			}*/
		return paretoRankSetListAll;
	}

	//拥挤距离计算
	public void crowdingDistance(ArrayList<ArrayList<ArrayList<Chromsome>>> paretoRankSetListAll) {
		//System.out.println("paretoRankSetListAll.size()="+paretoRankSetListAll.size());
		for(int i = 0;i < paretoRankSetListAll.size();i++) {
			for (int j = 0; j < paretoRankSetListAll.get(i).size(); j++) {
				ArrayList<Chromsome> paretoRankSet = new ArrayList<Chromsome>();
				paretoRankSet = paretoRankSetListAll.get(i).get(j);
				//System.out.println("paretoRankSet="+paretoRankSet.size());
				if (paretoRankSet.size() > 0 && paretoRankSet != null) {
					//按照目标依次计算个体的crowdingDistance值
					for (int m = 0; m < 3; m++) {
						//对paretoRankSet中的个体按照fitness从小到大进行排序
						paretoRankSet = this.getSortList(paretoRankSet, m);
						if (paretoRankSet.size() == 1) {// paretoRankSet中只有一个个体
							paretoRankSet.get(0).crowdingDistance += 100.0;
						} else if (paretoRankSet.size() == 2) {// paretoRankSet中有两个个体
							paretoRankSet.get(0).crowdingDistance += 100.0;
							paretoRankSet.get(paretoRankSet.size() - 1).crowdingDistance += 1000000.0;

						} else {// paretoRankSet个体数大于2
							double minFitness = paretoRankSet.get(0).fitness[m];// 目标值的最小值
							double maxFitness = paretoRankSet.get(paretoRankSet.size() - 1).fitness[m];// 目标值的最大值

							paretoRankSet.get(0).crowdingDistance += 100.0;// 将第1个和最后1个个体的Distance设置为无穷大，以便保存处在边界上的个体，为了便于计算设置为10000.0
							paretoRankSet.get(paretoRankSet.size() - 1).crowdingDistance += 1000000.0;

							for (int k = 1; k <= paretoRankSet.size() - 2; k++) {// 依次计算其他个体的crowdingDistance
								paretoRankSet.get(k).crowdingDistance += ((paretoRankSet.get(k + 1).fitness[m] - paretoRankSet.get(k - 1).fitness[m]) / (maxFitness - minFitness));
							}
						}
					}
				}
			}
		}
		//System.out.println();
		//System.out.println("paretoRankSetListAll.size()="+paretoRankSetListAll.size());
		// 在一个paretoRankSet中根据个体的crowdingDistance从大到小对个体进行排序
		for(int i = 0;i < paretoRankSetListAll.size();i++) {
			for (int j = 0; j < paretoRankSetListAll.get(i).size(); j++) {
				ArrayList<Chromsome> paretoRankSet = paretoRankSetListAll.get(i).get(j);
				//System.out.println("paretoRankSet="+paretoRankSet.size());
				// 对上面得到的paretoRankSet中的个体根据其crowdingDistance从大到小进行排序
				for (int m = 0; m < paretoRankSet.size(); m++) {
					Chromsome temp = new Chromsome(3);
					for (int k = m + 1; k < paretoRankSet.size(); k++) {
						if (paretoRankSet.get(m).crowdingDistance < paretoRankSet.get(k).crowdingDistance) {
							temp = paretoRankSet.get(m);
							paretoRankSet.set(m, paretoRankSet.get(k));
							paretoRankSet.set(k, temp);
						}
					}
				}
			}
		}
			/*System.out.println("排完序");
			for(int i = 0;i < init.getPopList().size();i++) {
				System.out.println("第" + i + "个：");
				for(int j = 0;j < paretoRankSetListAll.get(i).size();j++) {
					System.out.println("第" + j + "层");
					for(int m = 0;m < paretoRankSetListAll.get(i).get(j).size();m++) {
						for(int n = 0;n < b;n++) {
							System.out.print(paretoRankSetListAll.get(i).get(j).get(m).deviceId.get(n)+" ");
						}
						System.out.println();
						for(int l = 0;l < 3;l++) {
							System.out.print(paretoRankSetListAll.get(i).get(j).get(m).fitness[l]+" ");
						}
						System.out.print(paretoRankSetListAll.get(i).get(j).get(m).crowdingDistance);
						System.out.println();
					}
					System.out.println();
				}
				System.out.println();
			}*/
	}

	public void makeNewPopulation() {
		// 对两代种群快速非支配排序
		ArrayList<ArrayList<ArrayList<Chromsome>>> paretoRankSetListAll = this.fastNondominateSort();
		// 对每个pareto前沿计算拥挤距离并排序
		//System.out.println("OOOO");
		this.crowdingDistance(paretoRankSetListAll);
		// 记录加入到popList中的个体总数
		int numChromsome = 0;
		// popList的下标
		int popListIndex = 0;
		//取出前面n个个体作为下一代种群
		//System.out.println("BBBBBBBBBBBBB");
		//System.out.println("paretoRankSetListAll.size()="+paretoRankSetListAll.size());
		for(int i = 0;i < paretoRankSetListAll.size();i++) {
			for (int j = 0; j < paretoRankSetListAll.get(i).size(); j++) {
				ArrayList<Chromsome> paretoRankSet = paretoRankSetListAll.get(i).get(j);
				numChromsome += paretoRankSet.size();

				if (numChromsome < insidePopSize && popListIndex < insidePopSize) {
					for (int m = 0; m < paretoRankSet.size(); m++) {
						migrationPopList.get(i).set(popListIndex++, paretoRankSet.get(m));
						// popList.get(popListIndex++).setBinChrom(paretoRankSet.get(j).getBinChrom());
					}
				} else {
					int overNum = numChromsome - insidePopSize;
					for (int m = 0; m < paretoRankSet.size() - overNum && popListIndex < insidePopSize; m++) {
						migrationPopList.get(i).set(popListIndex++, paretoRankSet.get(m));

						// popList.get(popListIndex++).setBinChrom(paretoRankSet.get(j).getBinChrom());
					}
				}
			}
			numChromsome = 0;
			popListIndex = 0;
		}

			/*System.out.println("新一代");
			for(int i = 0;i < init.getPopList().size();i++) {
				for(int j = 0;j < 10;j++) {
					for(int l = 0;l < b;l++) {
						System.out.print(migrationPopList.get(i).get(j).deviceId.get(l)+" ");
					}
					System.out.println();
					for(int x = 0;x < 3;x++) {
						System.out.print(migrationPopList.get(i).get(j).fitness[x]+" ");
					}
					System.out.println();
				}
				System.out.println();
			}*/
	}

	public void runInsideNSGA2(){

		// 初始化种群
		//init.initPopulation();
		this.initInsidePopulation();

		this.evaluateInside();

		//计算当前和适应度
		this.calFitness();

		// 将种群保存起来供下次迭代
		this.copyPopListToPrePopList();

		int gen = 0;

		while (gen < 50) {

			//进入下一代
			gen += 1;
			//交叉变异操作
			this.crossoverOperator();
			this.mutationOperator();
			//System.out.println("内部第"+gen+"代");
			this.repair();
			// 计算目前背包重量和适应度
			this.calFitness();


			//找出支配高，拥挤距离大的个体形成新的种群
			//System.out.println("第"+gen+"代");
			this.makeNewPopulation();
			//System.out.println("--------------");
			//将种群保存起来供下次迭代
			this.copyPopListToPrePopList();


		}

			/*this.crossoverOperator();
			this.mutationOperator();

			// 计算目前背包重量和适应度
			this.calFitness();


			//找出支配高，拥挤距离大的个体形成新的种群
			System.out.println("第"+gen+"代");
			this.makeNewPopulation();

			//将种群保存起来供下次迭代
			this.copyPopListToPrePopList();
			choose();*/
			/*System.out.println("最后一代：");
			for(int i = 0;i < init.getPopList().size();i++) {
				for(int j = 0;j < 10;j++) {
					for(int l = 0;l < b;l++) {
						System.out.print(migrationPopList.get(i).get(j).deviceId.get(l)+" ");
					}
					System.out.println();
				}
				System.out.println();
			}*/

	}

	public ArrayList<Chromsome> choose() {
		ArrayList<Chromsome> finalPopList = new ArrayList<Chromsome>();
		double min_t = 0.0;
		double max_t = 0.0;
		double min_e = 0.0;
		double max_e = 0.0;
		double min_cbf = 0.0;
		double max_cbf = 0.0;
		for(int i = 0;i < migrationPopList.size();i++) {
			min_t = migrationPopList.get(i).get(0).fitness[0];
			max_t = migrationPopList.get(i).get(0).fitness[0];
			min_e = migrationPopList.get(i).get(0).fitness[1];
			max_e = migrationPopList.get(i).get(0).fitness[1];
			min_cbf = migrationPopList.get(i).get(0).fitness[2];
			max_cbf = migrationPopList.get(i).get(0).fitness[2];
			double[] u_t = new double[migrationPopList.get(i).size()];
			double[] u_e = new double[migrationPopList.get(i).size()];
			double[] u_cbf = new double[migrationPopList.get(i).size()];
			double[] u = new double[migrationPopList.get(i).size()];
			for(int j = 1;j < migrationPopList.get(i).size();j++) {
				if(migrationPopList.get(i).get(j).fitness[0] < min_t) {
					min_t = migrationPopList.get(i).get(j).fitness[0];
				}
				if(migrationPopList.get(i).get(j).fitness[0] > max_t) {
					max_t = migrationPopList.get(i).get(j).fitness[0];
				}
				if(migrationPopList.get(i).get(j).fitness[1] < min_e) {
					min_e = migrationPopList.get(i).get(j).fitness[1];
				}
				if(migrationPopList.get(i).get(j).fitness[1] > max_e) {
					max_e = migrationPopList.get(i).get(j).fitness[1];
				}
				if(migrationPopList.get(i).get(j).fitness[2] < min_cbf) {
					min_cbf = migrationPopList.get(i).get(j).fitness[2];
				}
				if(migrationPopList.get(i).get(j).fitness[2] > max_cbf) {
					max_cbf = migrationPopList.get(i).get(j).fitness[2];
				}
			}
			for(int m = 0;m < migrationPopList.get(i).size();m++) {
				if(min_t == max_t) {
					u_t[m] = 1;
				}
				else {
					u_t[m] = (max_t - migrationPopList.get(i).get(m).fitness[0]) / (max_t - min_t);
				}
				if(min_e == max_e) {
					u_e[m] = 1;
				}
				else {
					u_e[m] = (max_e - migrationPopList.get(i).get(m).fitness[1]) / (max_e - min_e);
				}
				if(min_cbf == max_cbf) {
					u_cbf[m] = 1;
				}
				else {
					u_cbf[m] = (max_cbf - migrationPopList.get(i).get(m).fitness[2]) / (max_cbf - min_cbf);
				}
				u[m] = 0.6 * u_t[m] + 0.3 * u_e[m] - 0.1 * u_cbf[m];
			}
			double max_u = u[0];
			int point = 0;
			for(int n = 1;n < migrationPopList.get(i).size();n++) {
				if(u[n] > max_u) {
					max_u = u[n];
					point = n;
				}
			}
			/*for(int x = 0;x < finalPopList.size();x++) {
				Chromsome chromsome = new Chromsome(3);
				finalPopList.add(chromsome);
			}*/
			finalPopList.add(migrationPopList.get(i).get(point));
		}
		/*System.out.println("最终");
		for(int i = 0;i < finalPopList.size();i++) {
			for(int j = 0;j < finalPopList.get(i).deviceId.size();j++) {
				System.out.print(finalPopList.get(i).deviceId.get(j)+" ");

			}
			//System.out.println(finalPopList.get(i).fitness[0]);
			System.out.println();
		}*/
		//System.out.println("-------");
		for(int i = 0;i < init.getPopList().size();i++) {
			for(int j = 0;j < finalPopList.get(i).deviceId.size();j++) {
				init.getPopList().get(i).mgDeviceId.set(j,finalPopList.get(i).deviceId.get(j));
			}
		}
		for(int i = 0;i < init.getPopList().size();i++) {
				init.getPopList().get(i).insidefitness[0] = finalPopList.get(i).fitness[0] ;
				init.getPopList().get(i).insidefitness[1] = finalPopList.get(i).fitness[1] ;
		}
		return finalPopList;
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

	public void calOutsideFitness() {
		ArrayList<Chromsome> OutsidePopList = choose();
		for(int i = 0;i < OutsidePopList.size();i++) {
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
			for(int j = 0;j < OutsidePopList.get(i).deviceId.size();j++) {
				for(int m = 0;m < OutsidePopList.get(i).deviceId.size();m++) {
					//temp.add(OutsidePopList.get(i).deviceId.get(m));
					if(OutsidePopList.get(i).deviceId.get(m).equals(OutsidePopList.get(i).deviceId.get(j))
							&& m != OutsidePopList.get(i).deviceId.size() - 1
							&& OutsidePopList.get(i).deviceId.get(m) != OutsidePopList.get(i).deviceId.get(m+1)
							&& !temp.contains(OutsidePopList.get(i).deviceId.get(j))) {
						e_trs += (0.2 * init.getService().get(init.getServiceReq().get(m)).getDt()
								/ calRij(init.getDevice().get(OutsidePopList.get(i).deviceId.get(m)), init.getDevice().get(OutsidePopList.get(i).deviceId.get(m + 1))));
						t_trs += init.getService().get(init.getServiceReq().get(m)).getDt()
								/ calRij(init.getDevice().get(OutsidePopList.get(i).deviceId.get(m)), init.getDevice().get(OutsidePopList.get(i).deviceId.get(m + 1)));
					}
					if(OutsidePopList.get(i).deviceId.get(m).equals(OutsidePopList.get(i).deviceId.get(j)) &&
							!temp.contains(OutsidePopList.get(i).deviceId.get(j))) {
						e_inv += init.getService().get(init.getServiceReq().get(m)).getInv();
						e_cmp += Math.pow(10, -26) * init.getDevice().get(OutsidePopList.get(i).deviceId.get(m)).getF()
								* init.getDevice().get(OutsidePopList.get(i).deviceId.get(m)).getF()
								* init.getService().get(init.getServiceReq().get(m)).getCr();
						t_cmp += init.getService().get(init.getServiceReq().get(m)).getCr()
								/ init.getDevice().get(OutsidePopList.get(i).deviceId.get(m)).getF();
					}
				}

				if(!temp.contains(OutsidePopList.get(i).deviceId.get(j))){
					temp.add(OutsidePopList.get(i).deviceId.get(j));
					e += (e_inv + e_cmp + e_trs);
					tempE.add((e_inv + e_cmp + e_trs));
					e_inv = 0.0;
					e_cmp = 0.0;
					e_trs = 0.0;
				}
			}
			for(int x = 0;x < temp.size();x++) {
				spt += evaluateLoc(init.getDevice().get(temp.get(x)), network.getIteRound());
				lbd += (init.getDevice().get(temp.get(x)).getRsd() - tempE.get(x)) / init.getDevice().get(temp.get(x)).getRsd();
			}
			t += t_cmp + t_trs;
			spt = spt / temp.size();
			lbn = lbd / temp.size();
			init.getPopList().get(i).fitness[0] = e;
			init.getPopList().get(i).fitness[1] = t;
			init.getPopList().get(i).fitness[2] = spt;
			init.getPopList().get(i).fitness[3] = lbn;
		}

	}

	public Map runOutside() {
		//init.initPopulation();
		Map nsga2Result = new HashMap();
		double max_u= 0.0;
		runInsideNSGA2();
		choose();
		calOutsideFitness();
		init.copyToPre();
		int gen = 0;
		while(gen < maxGen) {
			gen += 1;
			//System.out.println(gen);
			init.crossover();
			init.mutation();
			runInsideNSGA2();
			choose();
			calOutsideFitness();
			init.makeOutsideNewPopulation();
			init.copyToPre();
		}
		nsga2Result = init.chooseOutside();
		return nsga2Result;
	}

	public Map runNoMigrationNSGA2() {
		Map NoMigrationResult = new HashMap();

		double max_u = 0.0;
		// 计算当前和适应度
		init.calNoMigrationFitness();
		// 将种群保存起来供下次迭代
		// this.preInitInsidePopulation();
		init.copyToPre();
		int gen = 0;
		while (gen < 200) {
			// 进入下一代
			gen += 1;
			//System.out.println(gen);
			// 交叉变异操作
			init.crossover();
			init.mutation();
			// 计算目前背包重量和适应度
			init.calNoMigrationFitness();
			// 找出支配高，拥挤距离大的个体形成新的种群
			// System.out.println("第"+gen+"代");
			init.makeOutsideNewPopulation();
			// 将种群保存起来供下次迭代
			init.copyToPre();
		}
		//max_u = init.chooseOutside();
		//return max_u;
		NoMigrationResult = init.chooseOutside();
		return NoMigrationResult;
	}

	public static Map runAll(double dataIntensiveRatio){

		/*
		 * 任务数量50，设备数10000，服务数量150个。初始种群大小为50。
		 * 普通设备比例是百分之八十。
		 * 图十图十一横坐标是数据密集型服务的概率，概率从百分之二十开始，百分之八十结束。图十纵坐标是消耗能量（选出的最优个体的服务组合消耗的能量）。
		 * 图十一纵坐标是消耗时间（选出的最优个体的服务组合消耗的能量）。
		 */
		ArrayList<Double> disErgAndTime1 = new ArrayList<>();
		ArrayList<Double> disErgAndTime2 = new ArrayList<>();
		ArrayList<Double> disErgAndTime3 = new ArrayList<>();
		ArrayList<Double> disErgAndTime4 = new ArrayList<>();

		CreateNetwork network = new CreateNetwork(30, 5, dataIntensiveRatio, 0.8);//设备数量 服务数量 数据密集型服务比例 普通设备占比
		//CreateNetwork network = new CreateNetwork(10000, 150, i, 0.8);

		InitPopulation init = new InitPopulation(10, 5, network, network.getIteRound(),5000,10000);
		//InitPopulation init = new InitPopulation(50, 10, network, network.getIteRound());

		init.initPopulation(); // Output initial population
		Inside inside = new Inside(10,init,5,network,100,5000,10000);
		//Inside inside = new Inside(50,init,10,network,200);

		ArrayList<Chromsome> save = new ArrayList<Chromsome>();
		for(int j = 0;j < init.getPopList().size();j++) {
			Chromsome chromsome = new Chromsome(10);
			save.add(chromsome);
		}


		for(int j = 0;j < init.getPopList().size();j++) {
			save.set(j, (Chromsome)init.getPopList().get(j).clone());
		}
		Map nsga2Result = new HashMap();
		nsga2Result = inside.runOutside();
		System.out.println("-------------------------------------");
		System.out.println("NSGA2迁移");
		ArrayList<Chromsome> nsga2Solution = (ArrayList<Chromsome>) nsga2Result.get("bestSolution");
		Chromsome nsga2ChromsomeSolution = nsga2Solution.get(0);
		System.out.println(nsga2ChromsomeSolution.deviceId);
		System.out.println(nsga2ChromsomeSolution.mgDeviceId);
		System.out.println(nsga2Result.get("conErg"));
		System.out.println(nsga2Result.get("conTime"));
		System.out.println("-------------------------------------");


/*
		for(int j = 0;j < save.size();j++) {
			init.getPopList().set(j, save.get(j).clone());
		}
		Map noMigrationResult = new HashMap();
		noMigrationResult = inside.runNoMigrationNSGA2();
		System.out.println("-------------------------------------");
		System.out.println("无迁移");
		ArrayList<Chromsome> noMigrationSolution = (ArrayList<Chromsome>) noMigrationResult.get("bestSolution");
		Chromsome noMigrationChromsomeSolution = noMigrationSolution.get(0);
		System.out.println(noMigrationChromsomeSolution.deviceId);
		System.out.println(noMigrationChromsomeSolution.mgDeviceId);
		System.out.println(noMigrationResult.get("conErg"));
		System.out.println(noMigrationResult.get("conTime"));
		System.out.println("-------------------------------------");
*/



/*
		for(int j = 0;j < save.size();j++) {
			init.getPopList().set(j, save.get(j).clone());
		}
		RandomMig ran =  new RandomMig(10, 100, 0.65, 0.05, save, 5, init.getService(), init.getDevice(), init.getServiceReq(), network.getIteRound(),5000,10000);
		//RandomMig ran =  new RandomMig(50, 200, 0.65, 0.05, save, 10, init.getService(), init.getDevice(), init.getServiceReq(), network.getIteRound());
		Map randomResult = new HashMap();
		randomResult = ran.runNSGA2();
		System.out.println("-------------------------------------");
		System.out.println("随机");
		ArrayList<Chromsome> randomSolution = (ArrayList<Chromsome>) randomResult.get("bestSolution");
		Chromsome randomChromsomeSolution = randomSolution.get(0);
		System.out.println(randomChromsomeSolution.deviceId);
		System.out.println(randomChromsomeSolution.mgDeviceId);
		System.out.println(randomResult.get("conErg"));
		System.out.println(randomResult.get("conTime"));
		System.out.println("-------------------------------------");
*/



		for(int j = 0;j < save.size();j++) {
			init.getPopList().set(j, save.get(j).clone());
		}
		Map greedyResult = new HashMap();
		GreedyMig G1 = new GreedyMig(init,100,0.65,0.05,network,5000,10000);
		//GreedyMig G1 = new GreedyMig(init,200,0.65,0.05,network);
		greedyResult = G1.runNSGA2();
		System.out.println("-------------------------------------");
		System.out.println("贪婪");
		ArrayList<Chromsome> greedySolution = (ArrayList<Chromsome>) greedyResult.get("bestSolution");
		Chromsome greedyChromsomeSolution = greedySolution.get(0);
		System.out.println(greedyChromsomeSolution.deviceId);
		System.out.println(greedyChromsomeSolution.mgDeviceId);
		System.out.println(greedyResult.get("conErg"));
		System.out.println(greedyResult.get("conTime"));
		System.out.println("-------------------------------------");

		Map allResult = new HashMap();
/*		Map random = new HashMap();
		random.put("deviceId",randomChromsomeSolution.deviceId);
		random.put("conErg",randomResult.get("conErg"));
		random.put("conTime",randomResult.get("conTime"));*/
		Map greedy = new HashMap();
		greedy.put("deviceId",greedyChromsomeSolution.deviceId);
		greedy.put("conErg",greedyResult.get("conErg"));
		greedy.put("conTime",greedyResult.get("conTime"));
		Map nsga2 = new HashMap();
		nsga2.put("deviceId",nsga2ChromsomeSolution.deviceId);
		nsga2.put("conErg",nsga2Result.get("conErg"));
		nsga2.put("conTime",nsga2Result.get("conTime"));

		//allResult.put("random",random);
		allResult.put("greedy",greedy);
		allResult.put("nsga2",nsga2);

		return  allResult;

	}





	public static void main(String[] args) throws IOException{

		/*
		 * 任务数量50，设备数10000，服务数量150个。初始种群大小为50。
		 * 普通设备比例是百分之八十。
		 * 图十图十一横坐标是数据密集型服务的概率，概率从百分之二十开始，百分之八十结束。图十纵坐标是消耗能量（选出的最优个体的服务组合消耗的能量）。
		 * 图十一纵坐标是消耗时间（选出的最优个体的服务组合消耗的能量）。
		 */

		double dataIntensiveRatio = 0.6;
		System.out.println("dataIntensiveRatio = " + dataIntensiveRatio);

		CreateNetwork network = new CreateNetwork(30, 5, dataIntensiveRatio, 0.8);//设备数量 服务数量 数据密集型服务比例 普通设备占比
		//CreateNetwork network = new CreateNetwork(10000, 150, i, 0.8);

		InitPopulation init = new InitPopulation(10, 5, network, network.getIteRound(),5000,10000);
		//InitPopulation init = new InitPopulation(50, 10, network, network.getIteRound());

		init.initPopulation(); // Output initial population
		Inside inside = new Inside(10,init,5,network,100,5000,10000);
		//Inside inside = new Inside(50,init,10,network,200);

		ArrayList<Chromsome> save = new ArrayList<Chromsome>();
		for(int j = 0;j < init.getPopList().size();j++) {
			Chromsome chromsome = new Chromsome(10);
			save.add(chromsome);
		}

		ArrayList serviceFlag = InitPopulation.serviceFlag;


		for(int j = 0;j < init.getPopList().size();j++) {
			save.set(j, (Chromsome)init.getPopList().get(j).clone());
		}
		Map nsga2Result = new HashMap();
		nsga2Result = inside.runOutside();
		System.out.println("-------------------------------------");
		System.out.println("NSGA2迁移");
		ArrayList<Chromsome> nsga2Solution = (ArrayList<Chromsome>) nsga2Result.get("bestSolution");
		Chromsome nsga2ChromsomeSolution = nsga2Solution.get(0);
		System.out.println(nsga2ChromsomeSolution.deviceId);
		System.out.println(nsga2ChromsomeSolution.mgDeviceId);
		System.out.println(nsga2Result.get("conErg"));
		System.out.println(nsga2Result.get("conTime"));
		System.out.println("-------------------------------------");


		for(int j = 0;j < save.size();j++) {
			init.getPopList().set(j, save.get(j).clone());
		}
		Map noMigrationResult = new HashMap();
		noMigrationResult = inside.runNoMigrationNSGA2();
		System.out.println("-------------------------------------");
		System.out.println("无迁移");
		ArrayList<Chromsome> noMigrationSolution = (ArrayList<Chromsome>) noMigrationResult.get("bestSolution");
		Chromsome noMigrationChromsomeSolution = noMigrationSolution.get(0);
		System.out.println(noMigrationChromsomeSolution.deviceId);
		System.out.println(noMigrationChromsomeSolution.mgDeviceId);
		System.out.println(noMigrationResult.get("conErg"));
		System.out.println(noMigrationResult.get("conTime"));
		System.out.println("-------------------------------------");



/*
		for(int j = 0;j < save.size();j++) {
			init.getPopList().set(j, save.get(j).clone());
		}
		RandomMig ran =  new RandomMig(10, 100, 0.65, 0.05, save, 5, init.getService(), init.getDevice(), init.getServiceReq(), network.getIteRound(),5000,10000);
		//RandomMig ran =  new RandomMig(50, 200, 0.65, 0.05, save, 10, init.getService(), init.getDevice(), init.getServiceReq(), network.getIteRound());
		Map randomResult = new HashMap();
		randomResult = ran.runNSGA2();
		System.out.println("-------------------------------------");
		System.out.println("随机");
		ArrayList<Chromsome> randomSolution = (ArrayList<Chromsome>) randomResult.get("bestSolution");
		Chromsome randomChromsomeSolution = randomSolution.get(0);
		System.out.println(randomChromsomeSolution.deviceId);
		System.out.println(randomChromsomeSolution.mgDeviceId);
		System.out.println(randomResult.get("conErg"));
		System.out.println(randomResult.get("conTime"));
		System.out.println("-------------------------------------");
*/



		for(int j = 0;j < save.size();j++) {
			init.getPopList().set(j, save.get(j).clone());
		}
		Map greedyResult = new HashMap();
		GreedyMig G1 = new GreedyMig(init,100,0.65,0.05,network,5000,10000);
		//GreedyMig G1 = new GreedyMig(init,200,0.65,0.05,network);
		greedyResult = G1.runNSGA2();
		System.out.println("-------------------------------------");
		System.out.println("贪婪");
		ArrayList<Chromsome> greedySolution = (ArrayList<Chromsome>) greedyResult.get("bestSolution");
		Chromsome greedyChromsomeSolution = greedySolution.get(0);
		System.out.println(greedyChromsomeSolution.deviceId);
		System.out.println(greedyChromsomeSolution.mgDeviceId);
		System.out.println(greedyResult.get("conErg"));
		System.out.println(greedyResult.get("conTime"));
		System.out.println("-------------------------------------");

	}
}
