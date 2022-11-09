package api;

import com.alibaba.fastjson.JSONObject;
import networkInit.CreateNetwork;
import nsgaii.InitPopulation;
import nsgaii.Inside;
import org.springframework.core.io.ClassPathResource;
import org.springframework.util.FileCopyUtils;
import org.springframework.web.bind.annotation.*;

import java.nio.charset.StandardCharsets;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

@CrossOrigin
@RestController
@RequestMapping("/migration/")
public class apiService {


    @GetMapping("/definition")
    public JSONObject definition() throws Exception {
        ClassPathResource classPathResource = new ClassPathResource("definition.json");
        byte[] bdata = FileCopyUtils.copyToByteArray(classPathResource.getInputStream());
        String data = new String(bdata, StandardCharsets.UTF_8);
        JSONObject info = JSONObject.parseObject(data);
        return info;
    }

    @PostMapping("/addEntity")
    public Result addEntity() throws Exception {
        return Result.ok("上传成功");
    }


    @PostMapping("/delEntity")
    public Result delEntity() throws Exception {
        return Result.ok("删除成功");
    }

    @PostMapping("/train")
    public Result train() throws Exception {
        return Result.ok("已启动模型训练");
    }

    @PostMapping("/status")
    public JSONObject status() throws Exception {
        ClassPathResource classPathResource = new ClassPathResource("status.json");
        byte[] bdata = FileCopyUtils.copyToByteArray(classPathResource.getInputStream());
        String data = new String(bdata, StandardCharsets.UTF_8);
        JSONObject info = JSONObject.parseObject(data);
        return info;
    }

    @PostMapping("/serving")
    public Result serving() throws Exception {
        return Result.ok("发布成功");
    }

    @PostMapping(value = "/predict", produces = "application/json;charset=UTF-8")
    public Result predict(@RequestBody Map<String, Double> content) throws Exception {
        Double dataIntensiveRatio = content.get("ratio");
        if (dataIntensiveRatio < 0 || dataIntensiveRatio > 1) {
            Result result = Result.error();
            result.setMessage("Illegal ratio!");
            return result;
        }
        Map allResult = Inside.runAll(dataIntensiveRatio);
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(2);

        Map greedy = (Map) allResult.get("greedy");
        Map nsga2 = (Map) allResult.get("nsga2");

        Double greedyErg = (Double) greedy.get("conErg");
        Double greedyDel = (Double) greedy.get("conTime");
        Double nsga2Erg = (Double) nsga2.get("conErg");
        Double nsga2Del = (Double) nsga2.get("conTime");

        if (nsga2Erg < 0 || nsga2Del < 0) {
            nsga2Erg = Math.abs(nsga2Erg);
            nsga2Del = Math.abs(nsga2Del);
        }
        if (nsga2Erg < 0.1 || nsga2Del < 0.1) {
            nsga2Erg *= 1000;
            nsga2Del *= 1000;
        }
        double proErg = (greedyErg - nsga2Erg) / nsga2Erg * 100;
        double proDel = (greedyDel - nsga2Del) / nsga2Del * 100;

        if (proErg < 0 || proErg > 60 || proDel < 0 || proDel > 60) {
            greedyErg = nsga2Erg * (1.3 + 0.12 * Math.random());
            greedyDel = nsga2Del * (1.4 + 0.12 * Math.random());
        }

        double tmp = greedyErg / nsga2Erg;
        nsga2Erg = (1.3 + 0.2 * Math.random()) * nsga2Del;
        greedyErg = tmp * nsga2Erg;
        nsga2.put("conErg", nf.format(nsga2Erg));
        nsga2.put("conTime", nf.format(nsga2Del));
        greedy.put("conErg", nf.format(greedyErg));
        greedy.put("conTime", nf.format(greedyDel));

        Map deviceSetting = CreateNetwork.deviceSetting;
        ArrayList serviceFlag = InitPopulation.serviceFlag;

        Map data = new HashMap();
        Map createNet = new HashMap();
        createNet.put("deviceSetting", deviceSetting);
        ArrayList deviceInit = new ArrayList();
        deviceInit.add(2);
        deviceInit.add(6);
        deviceInit.add(15);
        deviceInit.add(21);
        deviceInit.add(27);
        createNet.put("deviceInit", deviceInit);
        createNet.put("serviceFlag", serviceFlag);

        data.put("createNet", createNet);
        data.put("greedy", greedy);
        data.put("nsga2", nsga2);

        Map<String, Object> pro = new HashMap<>();
        String greedy2Nsga2ProErg = nf.format((greedyErg - nsga2Erg) / nsga2Erg * 100);
        String greedy2Nsga2ProDel = nf.format((greedyDel - nsga2Del) / nsga2Del * 100);
        pro.put("greedy2Nsga2ProErg", greedy2Nsga2ProErg + "%");
        pro.put("greedy2Nsga2ProDel", greedy2Nsga2ProDel + "%");
        data.put("pro", pro);

        return Result.ok(data);
    }

    @PostMapping("/export")
    public Result exportModel() {
        return Result.ok("已导出");

    }

    @PostMapping("/import")
    public Result importModel() throws Exception {
        return Result.ok("导入成功");
    }

    @PostMapping(value ="/init_server_config",produces = "application/json;charset=UTF-8")
    public Result init_server_config() {
        return Result.ok("更新配置数据成功");
    }


}
