package nsgaii;

import java.util.ArrayList;




public class Chromsome implements Cloneable{
	public ArrayList<Integer> deviceId;//基因序列
	public double[] insidefitness;//内层适应度值数组
	public double[] fitness;//外层适应度值数组
	public int numDominated;//支配当前个体的个体数
	private int paretoRank;//个体的pareto等级
	public double crowdingDistance;//当前个体的拥挤距离
	//private ArrayList<Chromsome> dominatingList;//当前个体支配的个体的集合
	public ArrayList<Integer> mgDeviceId;//当前个体对应的迁移策略设备id集合

    public Chromsome clone() {   
    	Chromsome o = null;   
        try {   
            o = (Chromsome) super.clone();
            ArrayList<Integer> tempdeviceId = new ArrayList<Integer>();
            ArrayList<Integer> tempMgDeviceId = new ArrayList<Integer>();
            //ArrayList<Chromsome> tempdominatingList = new ArrayList<Chromsome>();
            for (int i = 0 ; i<this.deviceId.size();i++) {
            	tempdeviceId.add(this.deviceId.get(i));
            }
            for(int i = 0;i < this.mgDeviceId.size();i++) {
            	tempMgDeviceId.add(this.mgDeviceId.get(i));
            }

            o.deviceId = tempdeviceId;
            o.mgDeviceId = tempMgDeviceId;
            o.fitness = this.fitness.clone();
           // o.dominatingList = tempdominatingList;
        } catch (CloneNotSupportedException e) {   
            e.printStackTrace();   
        }   
        return o;   
    }   
	
	public Chromsome() {
		insidefitness = null;
		fitness = null;
		deviceId = null;
		numDominated=0;
		paretoRank=0;
		crowdingDistance=0.0;
		mgDeviceId = null;
	//	dominatingList=new ArrayList<Chromsome>();
	}

	public Chromsome(int numObjective) {//基因序列的长度，即用户选择的服务的个数；      目标函数的个数
		this.deviceId = new ArrayList<Integer>();
		this.fitness = new double[numObjective];
		this.insidefitness = new double[numObjective];
		this.numDominated = 0;
		this.paretoRank = 0;
		this.crowdingDistance = 0.0;
		this.mgDeviceId = new ArrayList<Integer>();
		//this.dominatingList = new ArrayList<Chromsome>();
	}

	/*public int[] getDeviceId() {
		return deviceId;
	}


	public void setDeviceId(int[] deviceId) {
		this.deviceId = deviceId;
	}*/

	/*public double[] getFitness() {
		return fitness;
	}


	public void setFitness(double[] fitness) {
		this.fitness = fitness;
	}*/

	/*public int getNumDominated() {
		return numDominated;
	}


	public void setNumDominated(int numDominated) {
		this.numDominated = numDominated;
	}*/

	public int getParetoRank() {
		return paretoRank;
	}

	public void setParetoRank(int paretoRank) {
		this.paretoRank = paretoRank;
	}

	/*public double getCrowdingDistance() {
		return crowdingDistance;
	}


	public void setCrowdingDistance(double crowdingDistance) {
		this.crowdingDistance = crowdingDistance;
	}*/

	/*	public ArrayList<Chromsome> getDominatingList() {
			return dominatingList;
		}


		public void setDominatingList(ArrayList<Chromsome> dominatingList) {
			this.dominatingList = dominatingList;
		}*/
}
