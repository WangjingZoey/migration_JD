package networkInit;

import java.util.ArrayList;

public class DeviceInit {
	private float x;//坐标
	private float y;
	private int id;
	private String name;
	private int flag;//设备类型
	private float spt;//空间限制
	private float rsd;//剩余能量
	private float f;//计算能力
	private float stg;//存储能力
	private float bnd;//带宽
	//public ArrayList<ServiceInit> service;
	private int N_Tsk;//设备容纳的计算任务的数量
	private int N_Cntr;//容纳的最大的容器数（依据计算能力）
	private int N_Cntr_Max;//容纳的最大容器数 5
	
	
	public DeviceInit(int id, String name) {
		this.id = id;
		this.name = name; 
		//this.flag = flag;
		//this.spt = spt;
		//this.rsd = rsd;
		//this.f = f;
		//this.stg = stg;
		//this.bnd = bnd;
		//this.service = service;
		//N_Tsk = n_Tsk;
		//N_Cntr = n_Cntr;
		//N_Cntr_Max = n_Cntr_Max;
	}
	
	public int getId() {
		return id;
	}
	public void setId(int id) {
		this.id = id;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}

	public float getX() {
		return x;
	}

	public void setX(float x) {
		this.x = x;
	}

	public float getY() {
		return y;
	}

	public void setY(float y) {
		this.y = y;
	}

	public int getFlag() {
		return flag;
	}

	public void setFlag(int flag) {
		this.flag = flag;
	}

	public float getSpt() {
		return spt;
	}

	public void setSpt(float spt) {
		this.spt = spt;
	}

	public float getRsd() {
		return rsd;
	}

	public void setRsd(float rsd) {
		this.rsd = rsd;
	}

	public float getF() {
		return f;
	}

	public void setF(float d) {
		this.f = d;
	}

	public float getStg() {
		return stg;
	}

	public void setStg(float stg) {
		this.stg = stg;
	}

	public float getBnd() {
		return bnd;
	}

	public void setBnd(float bnd) {
		this.bnd = bnd;
	}

	public int getN_Tsk() {
		return N_Tsk;
	}

	public void setN_Tsk(int n_Tsk) {
		N_Tsk = n_Tsk;
	}

	public int getN_Cntr() {
		return N_Cntr;
	}

	public void setN_Cntr(int n_Cntr) {
		N_Cntr = n_Cntr;
	}

	public int getN_Cntr_Max() {
		return N_Cntr_Max;
	}

	public void setN_Cntr_Max(int n_Cntr_Max) {
		N_Cntr_Max = n_Cntr_Max;
	}
	
	
}
