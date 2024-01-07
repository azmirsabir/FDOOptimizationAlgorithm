/**
 * @author Jaza Mahmood Abdullah  jazamahmood@gmail.com
 * Apache License Version 2.0, January 2004
 * http://www.apache.org/licenses/
 */


import java.util.ArrayList;

// this class hase Bee search agent properties and details
public class Bee {

    ArrayList xis = new ArrayList();
    ArrayList lastPast = new ArrayList();
    boolean isDominated = true;
    int number;
    Functions fitness ;

    public boolean isDominated() {
        return isDominated;
    }

    public void setIsDominated(boolean isDominated) {
        this.isDominated = isDominated;
    }

    public int getNumber() {
        return number;
    }

    public void setNumber(int number) {
        this.number = number;
    }

    public Bee(ArrayList p, int number, Functions fitness ) {
        this.xis.addAll(p);
        this.number = number;
        this.fitness =fitness;
    }

    public ArrayList getXis() {
        return xis;
    }

    public void setXis(ArrayList p) {
        this.xis = p;
    }

    public ArrayList getLastPace() {
        return lastPast;
    }

    public void setLastPast(ArrayList lastPast) {
        this.lastPast = lastPast;
    }

    @Override
    public String toString() {
        String str="";
        for(Object d : this.xis){
            str += (double)d+ ", ";
        }
        str = str.substring(0, str.length()-2);
        return str; //To change body of generated methods, choose Tools | Templates.
    }


    double getBeeFitness(String functionName) {
        double ft = 0.0;
        switch (functionName) {
            case "function_FM":
                ft = fitness.function_FM(xis);
                break;
            case "function_antenna":
                ft = fitness.function_antenna(xis);
                break;
            case "FX1":
                ft = fitness.function_1(xis);
                break;
            case "FX2":
                ft = fitness.function_2(xis);
                break;
            case "FX3":
                ft = fitness.function_3(xis);
                break;
            case "FX4":
                ft = fitness.function_4(xis);
                break;
            case "FX5":
                ft = fitness.function_5(xis);
                break;
            case "FX6":
                ft = fitness.function_6(xis);
                break;
            case "FX7":
                ft = fitness.function_7(xis);
                break;
            case "FX8":
                ft = fitness.function_8(xis);
                break;
            case "FX9":
                ft = fitness.function_9(xis);
                break;
            case "FX10":
                ft = fitness.function_10(xis);
                break;
            case "FX11":
                ft = fitness.function_11(xis);
                break;
            case "FX12":
                ft = fitness.function_12(xis);
                break;
            case "FX13":
                ft = fitness.function_13(xis);
                break;
            case "FX14":
                ft = fitness.function_14(xis);
                break;
            case "FX15":
                ft = fitness.function_15(xis);
                break;
            case "FX16":
                ft = fitness.function_16(xis);
                break;
            case "FX17":
                ft = fitness.function_17(xis);
                break;
            case "FX18":
                ft = fitness.function_18(xis);
                break;
            case "FX19":
                ft = fitness.function_19(xis);
                break;
            case "CEC01":
                ft = fitness.CEC01(xis);
                break;
            case "CEC02":
                ft = fitness.CEC02(xis);
                break;
            case "CEC03":
                ft = fitness.CEC03(xis);
                break;
            case "CEC04":
                ft = fitness.CEC04(xis);
                break;
            case "CEC05":
                ft = fitness.CEC05(xis);
                break;
            case "CEC06":
                ft = fitness.CEC06(xis);
                break;
            case "CEC07":
                ft = fitness.CEC07(xis);
                break;
            case "CEC08":
                ft = fitness.CEC08(xis);
                break;
            case "CEC09":
                ft = fitness.CEC09(xis);
                break;
            case "CEC10":
                ft = fitness.CEC10(xis);
                break;
        }
        return ft;

    }



}
