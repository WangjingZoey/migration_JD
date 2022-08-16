package networkInit;

import java.util.ArrayList;

public class ServiceInit {
	private int id;
	private String name;
	private int flag;//哪种服务 数据密集型 计算密集型
	private float wkd;//服务大小 未给 generateSkewnessService里初始化 1.0-2.0
	private int cr;//cpu占用量 500 1500
	private float dt;//内存占用量 
	private float bnd;//带宽 未给
	private int inv;//调用能耗  
	private ArrayList<DeviceInit> device;
	
	
	public ServiceInit(int id, ArrayList<DeviceInit> device) {
		
		this.id = id;
	//	this.name = name;
		this.device = device;
		//this.flag = flag;
		//this.wkd = wkd;
		//this.cr = cr;
		//this.stg = stg;
		//this.bnd = bnd;
		//this.device = device;
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
	public int getFlag() {
		return flag;
	}
	public void setFlag(int flag) {
		this.flag = flag;
	}
	public float getWkd() {
		return wkd;
	}
	public void setWkd(float d) {
		this.wkd = d;
	}
	public float getCr() {
		return cr;
	}
	public void setCr(int cr) {
		this.cr = cr;
	}
	
	public float getDt() {
		return dt;
	}

	public void setDt(float dt) {
		this.dt = dt;
	}

	public float getBnd() {
		return bnd;
	}
	public void setBnd(float bnd) {
		this.bnd = bnd;
	}
	
	public ArrayList<DeviceInit> getDevice() {
		return device;
	}

	public void setDevice(ArrayList<DeviceInit> device) {
		this.device = device;
	}

	public int getInv() {
		return inv;
	}

	public void setInv(int inv) {
		this.inv = inv;
	}
	
	
	
}
