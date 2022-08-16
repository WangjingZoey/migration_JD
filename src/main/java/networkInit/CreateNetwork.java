package networkInit;

import java.util.*;

import networkInit.Round;
import variable.Variable;
import sun.print.resources.serviceui;

public class CreateNetwork {
	private ArrayList<ServiceInit> service;
	private ArrayList<DeviceInit> device;
	private int numAllDevice;
	private int numAllService;
	double serviceSkewness;
	double deviceSkewness;
	Random random = new Random();
	Round iteRound;
	public static Map deviceSetting = new HashMap();
	//Random random = new Random();
	
	/*public CreateNetwork(ArrayList<ServiceInit> service, ArrayList<DeviceInit> device) {
		super();
		this.service = service;
		this.device = device;
	} */
	
	public ArrayList<ServiceInit> getService() {
		return service;
	}

	public CreateNetwork(int numAllDevice, int numAllService,double serviceSkewness,double deviceSkewness) {
		this.numAllDevice = numAllDevice;
		this.numAllService = numAllService;
		/*this.service = new ArrayList<ServiceInit>();
		this.device = new ArrayList<DeviceInit>();*/
		this.serviceSkewness = serviceSkewness;
		this.deviceSkewness = deviceSkewness;
		device = new ArrayList<DeviceInit>();
		service = new ArrayList<ServiceInit>();
	}

	public void setService(ArrayList<ServiceInit> service) {
		this.service = service;
	}

	public Round getIteRound() {
		return iteRound;
	}

	public void setIteRound(Round iteRound) {
		this.iteRound = iteRound;
	}

	public ArrayList<DeviceInit> getDevice() {
		return device;
	}

	public void setDevice(ArrayList<DeviceInit> device) {
		this.device = device;
	}

	public CreateNetwork() {
		this.service = new ArrayList<ServiceInit>();
		this.device = new ArrayList<DeviceInit>();
	} 

	public int getNumAllDevice() {
		return numAllDevice;
	}

	public void setNumAllDevice(int numAllDevice) {
		this.numAllDevice = numAllDevice;
	}

	public int getNumAllService() {
		return numAllService;
	}

	public void setNumAllService(int numAllService) {
		this.numAllService = numAllService;
	}

	private void preProcess() {
		Variable.cols = (int) (Variable.width / Variable.gSide);
		Variable.rows = (int) (Variable.height / Variable.gSide);
		if (Variable.height % Variable.gSide != 0) {
			Variable.rows += 1;
			Variable.height = (double) (Variable.rows * Variable.gSide);
		}
		if (Variable.width % Variable.gSide != 0) {
			Variable.cols += 1;
			Variable.width = (double) (Variable.cols * Variable.gSide);
		}
	}

	public void generateSkewnessDevice(ArrayList<DeviceInit> device) {
		preProcess();
		Random random = new Random(System.currentTimeMillis());
		float itex = (float) (Variable.width / 2 + Variable.startX);
		float itey = (float) (Variable.height / 2 + Variable.startY);
		float iter = Variable.height <= Variable.width ? (float) (Variable.height / 3) : (float) (Variable.width / 3);
		iteRound = new Round(itex, itey, iter);
		System.out.println(iteRound);
		// all
		for (int i = 0; i <  (int)(numAllDevice * (1 - Variable.skewness)); i++) {
			float x = (float) (random.nextFloat() * Variable.width + Variable.startX);
			float y = (float) (random.nextFloat() * Variable.height + Variable.startY);
			device.get(i).setX(x);
			device.get(i).setY(y);
			//device.get(i).setFlag(1);//普通设备
		}
		// left down corner
		int midX = Variable.cols / 2, midY = Variable.rows / 2;
		for (int i = (int) (numAllDevice * (1 - Variable.skewness)); i < ((int)  (numAllDevice * (1 - Variable.skewness)+numAllDevice * Variable.skewness / 2)); i++) {
			float x = (float) (random.nextFloat() * midX * Variable.gSide + Variable.startX);
			float y = (float) (random.nextFloat() * midY * Variable.gSide + Variable.startY);
			device.get(i).setX(x);
			device.get(i).setY(y);
			//device.get(i).setFlag(0);
		}
		// right upper
		for (int i = ((int)  (numAllDevice * (1 - Variable.skewness)+numAllDevice * Variable.skewness / 2)); i < numAllDevice; i++) {
			float x = (float) (random.nextFloat() * (Variable.cols - midX) * Variable.gSide + midX * Variable.gSide
					+ Variable.startX);
			float y = (float) (random.nextFloat() * (Variable.rows - midY) * Variable.gSide + midY * Variable.gSide
					+ Variable.startY);
			device.get(i).setX(x);
			device.get(i).setY(y);
			//device.get(i).setFlag(0);
			/*device.get(i).setRsd(Variable.initRsdErg);
			device.get(i).setStg(Variable.stg);
			device.get(i).setF((float) (random.nextFloat() * 800000 + 200000));*/
		}
		for(int i = 0;i < (int)(numAllDevice * (1 - deviceSkewness));i++){
			device.get(i).setFlag(0);
			device.get(i).setRsd(Variable.initRsdErg);							//剩余能量
			device.get(i).setStg(Variable.stg);									//存储能力
			device.get(i).setF((float) (random.nextFloat() * 800000 + 200000)); //计算能力200,000-1,000,000 //(float) (random.nextFloat() * 800000 + 200000)
		}
		for(int i = (int)(numAllDevice * (1 - deviceSkewness));i < numAllDevice;i++) {
			device.get(i).setFlag(1);//普通设备
			device.get(i).setRsd(Variable.initRsdErg);							//剩余能量
			device.get(i).setStg(Variable.stg);									//存储能力
			device.get(i).setF((float) (random.nextFloat() * 800000 + 200000)); //计算能力200,000-1,000,000  //(float) (random.nextFloat() * 800000 + 200000)
		}
		deviceSetting.clear();
		for(int i=0; i < numAllDevice; i++){
			Map currentDevice = new HashMap();
			currentDevice.put("id",device.get(i).getId());
			currentDevice.put("x",device.get(i).getX());
			currentDevice.put("y",device.get(i).getY());
			currentDevice.put("flag",device.get(i).getFlag());
			currentDevice.put("name",device.get(i).getName());
			deviceSetting.put(i,currentDevice);
		}

	}
	
	public void generateSkewnessService(ArrayList<ServiceInit> service) {
		int cnt0 = 0;
		if(serviceSkewness == 0) cnt0 = 0;
		else if(serviceSkewness > 0 && serviceSkewness <= 0.2) cnt0 = 1;
		else if(serviceSkewness > 0.2 && serviceSkewness <= 0.4) cnt0 = 2;
		else if(serviceSkewness > 0.4 && serviceSkewness <= 0.6) cnt0 = 3;
		else if(serviceSkewness > 0.6 && serviceSkewness <= 0.8) cnt0 = 4;
		else cnt0 = 5;
		int cnt1 = 5 - cnt0;

		for (int i = 0; i < cnt0; i++) {
			service.get(i).setFlag(0);//数据密集型
			service.get(i).setInv(1);
			service.get(i).setWkd((float) (random.nextFloat() + 1));
			service.get(i).setCr(random.nextInt(600) + 200);			//cpu占用量 200-800
			service.get(i).setDt(random.nextFloat() * 1024000  + 2048000 );	//内存占用量 2048000-3072000
		}
		for (int i = 0; i < cnt1; i++) {
			service.get(i).setFlag(1);//计算密集型
			service.get(i).setInv(2);
			service.get(i).setWkd((float) (random.nextFloat() + 1));
			service.get(i).setCr(random.nextInt(1000) + 1000);		//cpu占用量 1000-2000
			service.get(i).setDt(random.nextFloat() * 921600 + 102400);		//内存占用量 1024000-1945600
		}





/*		for (int i = 0; i <  (int)(numAllService * (1 - serviceSkewness)); i++) {
			service.get(i).setFlag(1);//计算密集型
			service.get(i).setInv(2);
			service.get(i).setWkd((float) (random.nextFloat() + 1));
			service.get(i).setCr(random.nextInt(1000) + 1000);		//cpu占用量 1000-2000
			service.get(i).setDt(random.nextFloat() * 921600 + 102400);		//内存占用量 1024000-1945600
		}
		for (int i = (int)(numAllService * (1 - serviceSkewness)); i < numAllService; i++) {
			service.get(i).setFlag(0);//数据密集型
			service.get(i).setInv(1);
			service.get(i).setWkd((float) (random.nextFloat() + 1));
			service.get(i).setCr(random.nextInt(600) + 200);			//cpu占用量 200-800
			service.get(i).setDt(random.nextFloat() * 1024000  + 2048000 );	//内存占用量 2048000-3072000
		}*/
	}

	public void initDevice(ArrayList<DeviceInit> device) {
		//device = new ArrayList<DeviceInit>();
		for(int i = 0;i < numAllDevice; i++) {
			DeviceInit deviceInit = new DeviceInit(i, i+"th device");
			device.add(deviceInit);
		}
	}
	
	public void initService(ArrayList<ServiceInit> service) {
		//service = new ArrayList<ServiceInit>();
		for(int i = 0;i < numAllService; i++) {
			ServiceInit serviceInit = new ServiceInit(i, new ArrayList<DeviceInit>());
			service.add(serviceInit);
		}
	}

	public void allocate(ArrayList<ServiceInit> service,ArrayList<DeviceInit> device) {
		int a = device.size();
		for(int i = 0;i < a;i++) {
			Random random = new Random();
			service.get(random.nextInt(service.size())).getDevice().add(device.get(i));
		}
	}
	
	
	
	
	/*public static void main(String[] args) {
		CreateNetwork createnetwork = new CreateNetwork();
		ArrayList<ServiceInit> service = new ArrayList<ServiceInit>();
		ArrayList<DeviceInit> device = new ArrayList<DeviceInit>();
		createnetwork.initDevice(device);
		createnetwork.initService(service);
		createnetwork.allocate(service,device);
		
		for(int i = 0;i < service.size();i++) {
			for(int j = 0;j < service.get(i).getDevice().size();j++) {
				System.out.print(service.get(i).getDevice().get(j).getId()+" ");
			}
			System.out.println();
		}
	}*/	
}
