package api;

import java.io.Serializable;

/**
 * @Description：后台接口返回json格式封装类
 */
public class Result<T> implements Serializable {

    //状态码
    private Integer code;
    //状态
    private String msg;
    //返回封装数据
    private T data;


    public Result() {
    }

    public Result(Integer code) {
        this.code = code;
    }

    public Result(Integer code, String msg) {
        this.code = code;
        this.msg = msg;
    }

    public Result(Integer code, String msg, T data) {
        this.code = code;
        this.msg = msg;
        this.data = data;
    }

    //不返回数据构造方法
    public Result(CodeEnum codeEnum) {
        this.code = codeEnum.getCode();
        this.msg = codeEnum.getMsg();
    }

    //返回数据构造方法
    public Result(CodeEnum codeEnum, T data) {
        this(codeEnum);
        this.data = data;
    }

    public Integer getCode() {
        return code;
    }

    public void setCode(Integer code) {
        this.code = code;
    }

    public String getMsg() {
        return msg;
    }

    public void setMsg(String msg) {
        this.msg = msg;
    }

    public T getData() {
        return data;
    }

    public void setData(T data) {
        this.data = data;
    }


    //请求成功（不返回数据）
    public static <T> Result <T> ok(){
        return new Result <T>(CodeEnum.SUCCESS);
    }

    //请求成功（返回数据）
    public static <T> Result <T> ok(T data){
        return new Result <T>(CodeEnum.SUCCESS,data);
    }

    //参数格式不正确
    public static <T> Result <T> error(){
        return new Result <T>(CodeEnum.BAD_REQUEST);
    }




    @Override
    public String toString() {
        return "Result{" +
                "code=" + code +
                ", msg='" + msg + '\'' +
                ", data=" + data +
                '}';
    }

}

