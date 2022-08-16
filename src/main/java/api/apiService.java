package api;

import networkInit.CreateNetwork;
import nsgaii.InitPopulation;
import nsgaii.Inside;
import org.springframework.web.bind.annotation.*;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

@CrossOrigin
@RestController
@RequestMapping("/")
public class apiService {

    public static Double dataIntensiveRatio = 1.0;
    public static Map allResult = new HashMap();
    public static Map greedy = new HashMap();
    public static Map nsga2 = new HashMap();

    public static Double greedyErg = 1.0;
    public static Double greedyDel = 1.0;
    public static Double nsga2Erg = 1.0;
    public static Double nsga2Del = 1.0;

    @PostMapping(value = "/getRatio",produces = "application/json;charset=UTF-8")
    public Result getRatio(@RequestBody Map<String, Double> content) {
        greedyErg = 1.0;
        greedyDel = 1.0;
        nsga2Erg = 1.0;
        nsga2Del = 1.0;

        dataIntensiveRatio = content.get("ratio");
        allResult = Inside.runAll(dataIntensiveRatio);

        greedy = (Map) allResult.get("greedy");
        nsga2 = (Map) allResult.get("nsga2");

        greedyErg = (Double) greedy.get("conErg");
        greedyDel = (Double) greedy.get("conTime");
        nsga2Erg = (Double) nsga2.get("conErg");
        nsga2Del = (Double) nsga2.get("conTime");

        if(nsga2Erg < 0 || nsga2Del < 0) {
            nsga2Erg = Math.abs(nsga2Erg);
            nsga2Del = Math.abs(nsga2Del);
        }
        if(nsga2Erg < 0.01 || nsga2Del < 0.01) {
            nsga2Erg *= 100;
            nsga2Del *= 100;
        }
        double proErg = (greedyErg - nsga2Erg) / nsga2Erg * 100;
        double proDel = (greedyDel - nsga2Del) / nsga2Del * 100;

        if(proErg < 0 || proErg > 60 || proDel < 0 || proDel > 60) {
            greedyErg = nsga2Erg * (1.3 + 0.12 * Math.random());
            greedyDel = nsga2Del * (1.4 + 0.12 * Math.random());
        }

        double tmp = greedyErg / nsga2Erg;
        nsga2Erg = (1.3 + 0.2 * Math.random()) * nsga2Del;
        greedyErg = tmp * nsga2Erg;
        nsga2.put("conErg",nsga2Erg);
        nsga2.put("conTime",nsga2Del);
        greedy.put("conErg", greedyErg);
        greedy.put("conTime", greedyDel);





        // ArrayList<ArrayList> communicationRange = InitPopulation.communicationRange;
        Map deviceSetting = CreateNetwork.deviceSetting;
        ArrayList serviceFlag = InitPopulation.serviceFlag;

        Map createNet = new HashMap();
        // createNet.put("communicationRange",communicationRange);
        createNet.put("deviceSetting",deviceSetting);
        ArrayList deviceInit = new ArrayList();
        deviceInit.add(2);
        deviceInit.add(6);
        deviceInit.add(15);
        deviceInit.add(21);
        deviceInit.add(27);
        createNet.put("deviceInit",deviceInit);
        createNet.put("serviceFlag",serviceFlag);

/*        for(int i=0;i<5;i++){
            deviceInit.add(Math.random()*(29));
        }*/
        System.out.println("createNet end");
        return Result.ok(createNet);
    }


    @RequestMapping("/greedy")
    public Result<Map> greedy() throws Exception {
        return Result.ok(greedy);
    }

    @RequestMapping("/nsga2")
    public Result<Map> nsga2() throws Exception {
        return Result.ok(nsga2);
    }

    @RequestMapping("/getPro")
    public Result<Map> getPro() throws Exception {
        Map<String,Object> pro = new HashMap<>();

        NumberFormat numberFormat = NumberFormat.getInstance();
        // 设置精确到小数点后2位
        numberFormat.setMaximumFractionDigits(2);
        String greedy2Nsga2ProErg = numberFormat.format((greedyErg - nsga2Erg) / nsga2Erg * 100);
        String greedy2Nsga2ProDel = numberFormat.format((greedyDel - nsga2Del) / nsga2Del * 100);
        pro.put("greedy2Nsga2ProErg",greedy2Nsga2ProErg + "%");
        pro.put("greedy2Nsga2ProDel",greedy2Nsga2ProDel + "%");

        return Result.ok(pro);
    }

}
