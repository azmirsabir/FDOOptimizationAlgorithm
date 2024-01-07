
/**
 * @author Jaza Mahmood Abdullah  jazamahmood@gmail.com
 * Apache License Version 2.0, January 2004
 * http://www.apache.org/licenses/
 */


import java.util.ArrayList;
import java.util.Random;


public class Functions {
    String functionName;

    public Functions(String functionName) {
        this.functionName = functionName;

    }

    public Double function_FM(ArrayList<Double> xs) {
        Double fitness = 0.0;
        Double y = 0.0;
        Double yn = 0.0;
        final double PI = Math.PI;
        double theta = (2 * PI) / 100;
        double a1 = xs.get(0);
        double w1 = xs.get(1);
        double a2 = xs.get(2);
        double w2 = xs.get(3);
        double a3 = xs.get(4);
        double w3 = xs.get(5);

        for (int i = 0; i < xs.size(); i++) {
            double x = xs.get(i).doubleValue();
            if (x > 6.3) {
                xs.set(i, 6.3 * (new Random().nextDouble()));
            } else if (x < -6.4) {
                xs.set(i, -6.4 * (new Random().nextDouble()));
            }
        }

        for (int t = 0; t <= 100; t++) {
            y = a1 * Math.sin(w1 * t * theta + a2 * Math.sin(w2 * t * theta + a3 * Math.sin(w3 * t * theta)));
            yn = 1.0 * Math.sin(5.0 * t * theta + 1.5 * Math.sin(4.8 * t * theta + 2.0 * Math.sin(4.9 * t * theta)));
            fitness += Math.pow(y - yn, 2);
        }
        return fitness;
    }

    public Double function_antenna(ArrayList<Double> xs) {
        Double fitness;
        Double sum = 0.0;
        final double PI = Math.PI;
        int evalationAngle = 45;
        int beamangle = 90;
        Random r = new Random();
        for (int j = 0; j < 4; j++) {
            double x = xs.get(j).doubleValue();
            if (x > 2.0) {
                xs.set(j, 2.0 * r.nextDouble());
            } else if (x <= 0.125) {
                xs.set(j, 0.125 + x * r.nextDouble());
            }
        }
        for (int i = 0; i < 4; i++) {
            double xi = xs.get(i).doubleValue();
            for (int j = i; j < 4; j++) {
                if (i != j) {
                    double xj = xs.get(j).doubleValue();
                    if (xi > xj && xi - xj < 0.26) {
                        xs.set(i, xi + 0.125);
                        xs.set(j, xj - 0.125);
                    } else if (xj > xi && xj - xi < 0.26) {
                        xs.set(j, xj + 0.125);
                        xs.set(i, xi - 0.125);
                    }
                }
            }
        }

        for (Double x : xs) {
            sum += Math.cos(2 * PI * x * (Math.cos(evalationAngle) - Math.cos(beamangle))) + Math.cos(2.25 * 2 * PI * (Math.cos(evalationAngle) - Math.cos(beamangle)));
        }
        fitness = 20 * Math.log(Math.abs(sum));
        return fitness;
    }

    public Double function_1(ArrayList<Double> xs) {
        Double fitness = 0.0;
        for (Double x : xs) {
            x = x - 30;
            fitness += Math.pow(x, 2);
        }

        return fitness;
    }

    public Double function_2(ArrayList<Double> xs) {
        Double fitness = 0.0;
        Double x_sum = 0.0;
        Double x_multiply = 1.0;
        for (double x : xs) {
            x = x - 3;
            x = Math.abs(x);
            x_multiply *= x;
            x_sum += x;
        }

        fitness = x_sum + x_multiply;
        return fitness;
    }

    public Double function_3(ArrayList<Double> xs) {
        Double fitness = 0.0;

        for (int i = 0; i < xs.size(); i++) {
            double sumXs = 0.0;
            for (int j = 0; j < i; j++) {
                sumXs += xs.get(j) - 30;
            }
            fitness += Math.pow(sumXs, 2);
        }

        return fitness;
    }

    public Double function_4(ArrayList<Double> xs) {
        Double fitness = Math.abs(xs.get(0));
        for (int i = 0; i < xs.size(); i++) {
            if ((double) Math.abs(xs.get(i)) > (double) fitness) {
                fitness = Math.abs(xs.get(i) - 30);
            }
        }
        return fitness;
    }

    public Double function_5(ArrayList<Double> xs) {
        Double fitness = 0.0;
        for (int i = 0; i < xs.size() - 1; i++) {
            double x = xs.get(i) - 15;
            fitness += (100 * Math.pow(((xs.get(i + 1) - 15) - Math.pow(x, 2)), 2)) + Math.pow((x - 1), 2);
        }
        return fitness;
    }

    public Double function_6(ArrayList<Double> xs) {
        Double fitness = 0.0;
        for (int i = 0; i < xs.size(); i++) {
            double x = xs.get(i) - 750;
            fitness += Math.pow(x + 0.5, 2);
        }
        return fitness;
    }

    public Double function_7(ArrayList<Double> xs) {
        Double fitness = 0.0;
        for (int i = 0; i < xs.size(); i++) {
            double x = xs.get(i) - 0.25;
            fitness += (i + 1) * Math.pow(x, 4);
        }
        return fitness + Math.random();
    }

    public Double function_8(ArrayList<Double> xs) {
        Double fitness = 0.0;
        for (int i = 0; i < xs.size(); i++) {
            double x = xs.get(i) - 300;
            fitness += Math.pow(-x, 2) * Math.sin(Math.sqrt(Math.abs(x)));
        }
        return fitness;
    }

    public Double function_9(ArrayList<Double> xs) {
        Double fitness = 0.0;
        for (int i = 0; i < xs.size(); i++) {
            double x = xs.get(i) - 2;
            fitness += Math.pow(x, 2) - (10 * Math.cos(2 * Math.PI * x)) + 10;
        }
        return fitness;
    }

    public Double function_10(ArrayList<Double> xs) {
        Double fitness = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int i = 0; i < xs.size(); i++) {
            sum1 += Math.pow(xs.get(i), 2);
            sum2 += Math.cos(2 * Math.PI * xs.get(i));
        }
        double n = xs.size();
        sum1 = -0.2 * Math.sqrt(sum1 / n);
        fitness = -20 * Math.exp(sum1) - Math.exp(sum2 / n) + 20 + Math.E;
        return fitness;
    }

    public Double function_11(ArrayList<Double> xs) {
        Double fitness = 0.0;
        double sum = 0.0;
        double prod = 1;
        for (int i = 0; i < xs.size(); i++) {
            double x = xs.get(i) - 400;
            sum += Math.pow(x, 2);
            prod *= Math.cos(x / Math.sqrt(i + 1));
            //    
        }
        //   System.out.print("prod" +prod);
        fitness = (1 / 4000) * sum - prod + 1;
        return fitness;
    }

    public Double function_12(ArrayList<Double> xs) {
        Double fitness = 0.0;
        double sum = 0.0;
        double u = 0.0;
        double yn = 0.0;
        for (int i = 0; i < xs.size(); i++) {
            double x = xs.get(i) - 30;
            double y = 1 + ((x + 1) / 4);
            sum += Math.pow((y - 1), 2) * Math.pow((1 + (10 * Math.sin(Math.PI * y + 1))), 2);
            u += Ufunc(x, 10, 100, 4);
        }
        yn = 1 + (((xs.get(xs.size() - 1) - 30) + 1) / 4);
        fitness = (Math.PI / xs.size()) * ((10 * Math.pow(Math.sin(Math.PI * (1 + ((xs.get(0) - 30) + 1) / 4)), 2)) + sum + Math.pow(yn - 1, 2)) + u;
        return fitness;
    }

    public Double function_13(ArrayList<Double> xs) {
        Double fitness = 0.0;
        double sum = 0.0;
        double sum2 = 0.0;
        double x_1 = 100 - xs.get(0);
        double x_n = 100 - xs.get(xs.size() - 1);
        final double PI = Math.PI;
        double u = 0.0;
        for (int i = 0; i < xs.size(); i++) {
            double x = xs.get(i) - 100;
            sum2 = Math.pow(Math.sin(3 * PI * x + 1), 2);
            sum += Math.pow(x - 1, 2) * (1 + sum2);
            u += Ufunc(x, 5, 100, 4);
        }
        fitness = 0.1 * ((Math.pow(Math.sin(3 * PI * (x_1)), 2)) +
                sum + Math.pow(x_n - 1, 2) *
                (1 + Math.pow(Math.sin(2 * PI * (x_n)), 2))) + u;
        return fitness;
    }

    private double Ufunc(double x, double a, double k, double m) {
        double u = 0.0;
        if (x > a) {
            u = k * Math.pow((x - a), m);
        } else if (x < a && x > -a) {
            u = 0.0;
        } else if (x < -a) {
            u = k * Math.pow((-x - a), m);
        }

        return u;
    }

    public Double function_14(ArrayList<Double> xs) {
        Double fitness = 0.0;

        for (int i = 0; i < xs.size(); i++) {
            double x = (1 + xs.get(i)) * 0.05;
            fitness += Math.pow(x, 2);
        }

        return fitness;
    }

    public Double function_15(ArrayList<Double> xs) {
        Double fitness = 0.0;
        double sum = 0.0;
        double mul = 1.0;

        for (int i = 0; i < xs.size(); i++) {
            double x = (1 + xs.get(i)) * 0.05;
            sum += Math.pow(x, 2) / 400;
            mul *= Math.cos(x / Math.sqrt(Math.abs(x)));
        }
        fitness = sum - mul + 1;
        return fitness;
    }

    public Double function_16(ArrayList<Double> xs) {
        Double fitness = 0.0;
        double sum = 0.0;
        double mul = 1.0;
////        double x1 = xs.get(0);
////        double x2 = xs.get(1);
////        fitness = 4 * Math.pow(x1, 2) - 2.1 * Math.pow(x1, 4) +(Math.pow(x1, 6) / 3) + x1 * x2 - 4 * Math.pow(x2, 2) + 4 * Math.pow(x2, 4);
//        double sum = 0.0;
//        double bee.multi = 1.0;
//        for (int i = 0; i < xs.size(); i++) {
//            double x = xs.get(i);
//            sum += Math.pow(x, 2) / 4000;
//            bee.multi *= Math.cos(x / Math.sqrt(i + 1));
//        }
//        fitness = sum - bee.multi + 1;
//        return fitness;
//    }
//
//    public Double function_17(ArrayList<Double> xs) {
//        Double fitness = 0.0;
//        double x1 = xs.get(0);
//        double x2 = xs.get(1);
//        fitness = Math.pow((x2 - (Math.pow(x1, 2)) * 5.1 / (4 * Math.pow(Math.PI, 2)) + 5 / Math.PI * x1 - 6), 2) + 10 * (1 - 1 / (8 * Math.PI)) * Math.cos(x1) + 10;
        for (int i = 0; i < xs.size(); i++) {
            double x = (1 + xs.get(i));
            sum += Math.pow(x, 2) / 400;
            mul *= Math.cos(x / Math.sqrt(Math.abs(x)));
        }
        fitness = sum - mul + 1;
        return fitness;
    }

    public Double function_17(ArrayList<Double> xs) {
        Double fitness = 0.0;
        double x1 = ((1 + xs.get(0)) * (5 / 32));
        double x2 = ((1 + xs.get(1)) * (5 / 32));
        double x3 = (1 + xs.get(2));
        double x4 = (1 + xs.get(3));
        double x5 = (1 + xs.get(4) * (5 / 0.5));
        double x6 = (1 + xs.get(5) * (5 / 0.5));
        double x7 = (1 + xs.get(6) * (5 / 100));
        double x8 = (1 + xs.get(7) * (5 / 100));
        double x9 = (1 + xs.get(8) * (5 / 100));
        double x10 = (1 + xs.get(9) * (5 / 100));

        double a = 20;
        double b = 0.2;
        double c = 2 * Math.PI;
        double sum = Math.pow(x1, 2) + Math.pow(x2, 2);

        double sumCos = Math.cos(Math.abs(x1) * c) + Math.cos(Math.abs(x2) * c);
        //  System.out.println("sumCos "+sumCos);
        double ackley = -a * Math.exp(-b * Math.sqrt(Math.abs(1 / 2 * sum))) - Math.exp(1 / 2 * sumCos) + a + Math.exp(1);

        sum = Math.pow(x3 - 10 * Math.cos(2 * Math.PI * x3), 2) + Math.pow(x4 - 10 * Math.cos(2 * Math.PI * x4), 2);
        double rastrigin = 20 + sum;

        double weierstrass = 0.5 * Math.sin(2 * x5) + 0.25 * Math.sin(4 * x6);

        sum = Math.pow(x7, 2) / 400 + Math.pow(x8, 2) / 400;
        double mul = Math.cos(x7 / Math.sqrt(Math.abs(x7))) * Math.cos(x8 / Math.sqrt(Math.abs(x8)));
        double griewank = sum - mul + 1;

        double sphere = Math.pow(x9, 2) + Math.pow(x10, 2);

        fitness = ackley + rastrigin + weierstrass + griewank + sphere;
        return fitness;

    }

    public Double function_18(ArrayList<Double> xs) {
        Double fitness = 0.0;
        double x1 = ((1 + xs.get(0)) * (1 / 5));
        double x2 = ((1 + xs.get(1)) * (1 / 5));
        double x3 = (1 + xs.get(2) * 10);
        double x4 = (1 + xs.get(3) * 10);
        double x5 = (1 + xs.get(4) * (0.05));
        double x6 = (1 + xs.get(5) * (0.05));
        double x7 = (1 + xs.get(6) * (5 / 32));
        double x8 = (1 + xs.get(7) * (5 / 100));
        double x9 = (1 + xs.get(8) * (5 / 100));
        double x10 = (1 + xs.get(9) * (5 / 100));
        double sum = 0.0;


        sum = Math.pow(x1 - 10 * Math.cos(2 * Math.PI * x1), 2) + Math.pow(x2 - 10 * Math.cos(2 * Math.PI * x2), 2);
        double rastrigin = 20 + sum;

        double weierstrass = 0.5 * Math.sin(2 * x3) + 0.25 * Math.sin(4 * x4);

        sum = Math.pow(x5, 2) / 400 + Math.pow(x5, 2) / 400;
        double mul = Math.cos(x5 / Math.sqrt(Math.abs(x5))) * Math.cos(x6 / Math.sqrt(Math.abs(x6)));
        double griewank = sum - mul + 1;

        double a = 20;
        double b = 0.2;
        double c = 2 * Math.PI;
        sum = Math.pow(x7, 2) + Math.pow(x8, 2);
        double sumCos = Math.cos(Math.abs(x7) * c) + Math.cos(Math.abs(x8) * c);
        double ackley = -a * Math.exp(-b * Math.sqrt(Math.abs(1 / 2 * sum))) - Math.exp(1 / 2 * sumCos) + a + Math.exp(1);


        double sphere = Math.pow(x9, 2) + Math.pow(x10, 2);

        fitness = ackley + rastrigin + weierstrass + griewank + sphere;
        return fitness;
    }

    public Double function_19(ArrayList<Double> xs) {
        Double fitness = 0.0;
//        double[][] aH = {
//            {3, 10, 30},
//            {0.1, 10, 35},
//            {3, 10, 30},
//            {0.1, 10, 35}
//        };
//        double[] cH = {1, 1.2, 3, 3.2};
//        double[][] pH = {
//            {0.3689, 0.117, 0.2673},
//            {0.4699, 0.4387, 0.747},
//            {0.1091, 0.8732, 0.5547},
//            {0.03815, 0.5743, 0.8828}
//        };
//
//        for (int i = 0; i < 4; i++) {
//            double sum = 0.0;
//            for (int j = 0; j < 3; j++) {
//                double x = xs.get(j);
//                sum += aH[i][j] * Math.pow((x - pH[i][j]), 2);
//            }
//            fitness = fitness - cH[i] * Math.exp(-sum);
//        }

        double x1 = (0.1 + xs.get(0) * (0.1 * 1 / 5));
        double x2 = (0.2 + xs.get(1) * (0.2 * 1 / 5));
        double x3 = (0.3 + xs.get(2) * 0.3 * 10);
        double x4 = (0.4 + xs.get(3) * 0.4 * 10);
        double x5 = (0.5 + xs.get(4) * (0.5 * 0.05));
        double x6 = (0.6 + xs.get(5) * (0.6 * 0.05));
        double x7 = (0.7 + xs.get(6) * 0.7 * (5 / 32));
        double x8 = (0.8 + xs.get(7) * 0.8 * (5 / 100));
        double x9 = (0.9 + xs.get(8) * 0.9 * (5 / 100));
        double x10 = (1 + xs.get(9) * (5 / 100));
        double sum = 0.0;


        sum = Math.pow(x1 - 10 * Math.cos(2 * Math.PI * x1), 2) + Math.pow(x2 - 10 * Math.cos(2 * Math.PI * x2), 2);
        double rastrigin = 20 + sum;

        double weierstrass = 0.5 * Math.sin(2 * x3) + 0.25 * Math.sin(4 * x4);

        sum = Math.pow(x5, 2) / 400 + Math.pow(x5, 2) / 400;
        double mul = Math.cos(x5 / Math.sqrt(Math.abs(x5))) * Math.cos(x6 / Math.sqrt(Math.abs(x6)));
        double griewank = sum - mul + 1;

        double a = 20;
        double b = 0.2;
        double c = 2 * Math.PI;
        sum = Math.pow(x7, 2) + Math.pow(x8, 2);
        double sumCos = Math.cos(Math.abs(x7) * c) + Math.cos(Math.abs(x8) * c);
        double ackley = -a * Math.exp(-b * Math.sqrt(Math.abs(1 / 2 * sum))) - Math.exp(1 / 2 * sumCos) + a + Math.exp(1);


        double sphere = Math.pow(x9, 2) + Math.pow(x10, 2);

        fitness = ackley + rastrigin + weierstrass + griewank + sphere;

        return fitness;
    }

    public double[] function_20(ArrayList<Double> xs) {
        double[] fx = {0.0, 0.0};
        int n = xs.size();
        fx[0] = xs.get(0);
        double sum = 0.0;
        for (int i = 1; i < n; i++) {
            sum += xs.get(i);
        }
        double gx = 1 + 9 / (n - 1) * sum;
        double h = 1 - Math.sqrt(Math.abs(fx[0] / gx));
        fx[1] = gx * h;
        return fx;
    }

    public double[] function_21(ArrayList<Double> xs) {
        double[] fx = {0.0, 0.0};
        int n = xs.size();
        fx[0] = xs.get(0);
        double sum = 0.0;
        for (int i = 1; i < n; i++) {
            sum += xs.get(i);
        }
        double gx = 1 + 9 / (n - 1) * sum;
        double h = 1 - Math.pow(fx[0] / gx, 2);
        fx[1] = gx * h;
        return fx;
    }

    public double[] function_22(ArrayList<Double> xs) {
        double[] fx = {0.0, 0.0};
        int n = xs.size();
        fx[0] = xs.get(0);
        double sum = 0.0;
        for (int i = 1; i < n; i++) {
            sum += xs.get(i);
        }
        double gx = 1 + (9 / 29) * sum;
        double h = 1 - Math.sqrt(Math.abs(fx[0] / gx)) - (fx[0] / gx) * Math.sin(10 * Math.PI * fx[0]);
        fx[1] = gx * h;
        return fx;
    }

    public double[] function_23(ArrayList<Double> xs) {
        double[] fx = {0.0, 0.0};
        int n = xs.size();
        fx[0] = xs.get(0);
        double sum = 0.0;
        for (int i = 1; i < n; i++) {
            sum += xs.get(i);
        }
        double gx = 1 + (9 / (n - 1)) * sum;
        double h = 1 - fx[0] / gx;
        fx[1] = gx * h;
        return fx;
    }

    public double[] function_24(ArrayList<Double> xs) {
        double[] fx = {0.0, 0.0, 0.0};
        int n = xs.size();
        fx[0] = xs.get(0);
        fx[1] = xs.get(1);
        double sum = 0.0;
        for (int i = 1; i < n; i++) {
            sum += xs.get(i);
        }
        double gx = 1 + (9 / (n - 1)) * sum;
        double h1 = 1 - Math.pow(fx[0] / gx, 2);
        double h2 = 1 - Math.pow(fx[1] / gx, 2);
        fx[2] = gx * h1 * h2;
        return fx;
    }

    public double CEC01(ArrayList<Double> xs) {
        double p1 = 0.0;
        double p2 = 0.0;
        double p3 = 0.0;
        double d = 72.661;
        double u = 0.0;
        double v = 0.0;
        double wk = 0.0;
        double pk = 0.0;
        int dim = xs.size();
        int m = (32 * dim);

        for (int j = 1; j <= dim; j++) {
            double x = xs.get(j - 1);
            u += x * Math.pow(1.2, dim - j);
        }
        if (u < d) u = Math.pow(u - d, 2);
        else u = 0;
        p1 = u;

        for (int j = 1; j <= dim; j++) {
            double x = xs.get(j - 1);
            v += x * Math.pow(-1.2, dim - j);
        }

        if (v < d) v = Math.pow(v - d, 2);
        else v = 0;
        p2 = v;


        for (int k = 0; k <= m; k++) {
            wk = 0.0;
            for (int j = 1; j <= xs.size(); j++) {
                double x = xs.get(j - 1);
                wk += x * Math.pow((2 * k / m) - 1, xs.size() - j);
            }
            if (wk > 1) pk += Math.pow(wk - 1, 2);
            else if (wk < 1) pk += Math.pow(wk + 1, 2);
            else pk += 0;
        }
        p3 = pk;

        return p1 + p2 + p3+ 1;
    }

    public double CEC02(ArrayList<Double> xs){
        int D = xs.size();
        int n = (int)Math.sqrt(D);
        double W = 0.0;
        double I =0;
        double H = 0.0;
        double Z = 0.0;

        for(int i=1; i <=n; i++){
            double x = xs.get(i-1);
            for(int k=1; k <= n; k++){
                if(i == k) I =1;
                else I =0;
                H = 1/(i+k-1);
                Z = x+(n*(k-1));
                W += Math.abs(H*Z-I);
            }
        }
        return W +1;
    }

    public double CEC03(ArrayList<Double> xs){
        int n = xs.size()/3;
        double d = 0.0;
        double sum =0.0;
        for(int i=1; i <=n-1; i++){
            double xi = xs.get(3*i-1);
            for(int j=i+1; j<=n; j++){
                double tmp=0.0;
                double xj = xs.get(3*j-1);
                for(int k=0; k <=2; k++){
                    tmp += Math.pow(xi+k-2-xj+k-2, 2);
                }
                d += Math.pow(tmp, 3);
                sum += (1/ Math.pow(d,2))- (2/d);
            }
        }
        return 12.7120622568+sum +1;
    }

    public double CEC04(ArrayList<Double> xs){
        int dim = xs.size();
        double x[][] = new double[1][dim];
        double sum =0.0;
        double shiftedMatrix[] = {4.3453613502650342e+01,  -7.5117860955706732e+01,   5.4110917436941946e+01,   2.1893626834216349e+00,  -3.3813797325740467e+00,  -3.0849165372014589e+01,   7.8077592550813023e+01,  -6.9901998485392895e+01,   3.7111456001695004e+01,   5.2241020487733664e+01};
        double rotatedMatrix[][] = {
                {8.8970810825119684e-01, 1.9871231543356224e-01, 3.5531377300377703e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,  -2.0660353462835387e-01, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {1.0419879983757413e-01,  -6.6358499459221376e-01, 4.5164451523757104e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 5.8720932972857365e-01, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {-4.3941933258454113e-01, 3.4165627723133662e-01, 8.1471256710105333e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,  -1.6255790164428213e-01, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 4.2863677588379034e-01, 0.0000000000000000e+00,  -6.3119168084438271e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 6.4642677573936591e-01, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 6.8398421839009127e-01, 0.0000000000000000e+00,  -1.0280182935971671e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 9.1093354741887778e-01},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.4071437781392207e-01, 0.0000000000000000e+00, 7.5339632462828454e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 6.4233436924473619e-01, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,  -5.2127926340507358e-01, 0.0000000000000000e+00, 1.8558504533667268e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,  -1.2500979819172911e-01},
                {6.6878564394300288e-02, 6.3516876408991430e-01,  -7.7542166446725458e-02, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 7.6555542658372966e-01, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,  -8.9245166717105195e-01, 0.0000000000000000e+00, -1.8436659152198701e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 4.1175111620270627e-01, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,  -1.2215803943548182e+00, 0.0000000000000000e+00,  -1.0421383218771203e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 4.5968278525834394e-01}
        };
        //shifting
        for (int j = 0; j < dim; j++) {
            x[0][j] = xs.get(j) - shiftedMatrix[j];
        }
        //rotating
        double shitedRotaatedX[][] = multiplyByMatrix(x, rotatedMatrix);

        for(int i =0; i < dim; i++ ){
            double Xi = shitedRotaatedX[0][i];
            sum += Math.pow(Xi, 2)- 10* Math.cos(2*Math.PI*Xi)+10;
        }
        return sum +1;
    }

    public double CEC05(ArrayList<Double> xs) {
        double sum = 0.0;
        double multi = 1.0;
        int i = 1;
        double x[][] = new double[1][xs.size()];
        double shiftedMatrix[] = {-1.6799910337105352e+01, 4.3906964270354706e+01, 2.4348491851402670e+01, -5.4897453475230122e+01, 5.8499441807390866e+01, 1.1845681821854726e-01, 7.0903743799265357e+01, -7.7796574718223610e-01, 4.4729687108066713e+01, -6.8148774722660320e+01};
        double rotatedMatrix[][] = {
                {-7.5988949123997229e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 5.9790648917707112e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -2.5509957135010197e-01, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {0.0000000000000000e+00, -6.4335422234689021e-02, 1.3912090644901491e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -7.0354948113015447e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 1.0874018967981698e+00, -9.3628657345778921e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -2.2957124131927584e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {-4.8529941749281533e-02, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.3915506135367807e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 9.3947788111907859e-01, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -7.3880347589611128e-01, -1.8005470643130739e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, -7.5105330032914350e-02, 6.4506504795749642e-01},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 2.2733754217810664e-01, 1.4610614422407869e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 8.7440899428783569e-01, 4.0296345646339404e-01},
                {0.0000000000000000e+00, -9.1794081957250695e-01, 1.0235511548555627e-02, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -9.3534397252777979e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {6.4823823233196287e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 7.2627933645267095e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -2.2870399993222773e-01, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -9.1460850041609643e-02, -8.8922500528071025e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.1795773407187566e-01, -3.1593778222470387e-01},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -6.2779134975715278e-01, 3.9435033573642242e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.5870737301184835e-01, -5.6718150042257454e-01}
        };
        //shifting
        for (int j = 0; j < xs.size(); j++) {
            x[0][j] = xs.get(j) - shiftedMatrix[j];
        }
        //rotating
        double shitedRotaatedX[][] = multiplyByMatrix(x, rotatedMatrix);

        for (int m = 0; m < xs.size(); m++) {
            sum += Math.pow(shitedRotaatedX[0][m], 2) / 4000;
            multi *= Math.cos(shitedRotaatedX[0][m] / Math.sqrt(i));
            i++;
        }
        return (sum - multi + 1) + 1;
    }

    public double CEC06(ArrayList<Double> xs){
        double a =0.5;
        double b =3.0;
        int kMax = 20;
        int D = xs.size();
        double sum = 0.0;
        double x[][] = new double[1][D];
        double shiftedMatrix[] = {4.4867071194977996e+01, 8.6557399521842626e-01,  -1.2297862364117918e+01,   2.9827246270062048e+01,   2.6528060932889602e+01,  -6.2879900924339843e+01,  -2.2494835379763892e+01,   9.3017723082107295e+00,   1.4887184097844738e+01,  -3.1096867523666873e+01 };
        double rotatedMatrix[][] = {
                {-1.5433743057196678e-01, 0.0000000000000000e+00, 7.7666311726871273e-01, 0.0000000000000000e+00, 1.1571979400226866e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00},
                { 0.0000000000000000e+00, 4.6806840267259536e-02, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,  -5.9264454599472804e-01, 1.6314935476659614e-01, 7.8737783169590370e-01, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {-1.7410812843826278e+00, 0.0000000000000000e+00,  -4.4194799352318298e-01, 0.0000000000000000e+00, 4.4605580480878959e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 2.7077411154472419e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,  -8.8999649318267127e-01, 3.6686185770629254e-01},
                {5.4888525059737507e-02, 0.0000000000000000e+00, 1.5570674387300532e+00, 0.0000000000000000e+00,  -3.0216546520289828e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 6.1164921138202333e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,  -6.1748299284504526e-01,  -1.5999277506278717e-01,  -4.6797682388189477e-01, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {0.0000000000000000e+00, -1.1226733726002835e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,  -1.3517591002752971e-01, 9.4075663040175728e-01,  -2.9000082877106131e-01, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 7.8172271740335475e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 4.9921405128267116e-01, 2.5052257846765580e-01, 2.7736863877405393e-01, 0.0000000000000000e+00, 0.0000000000000000e+00},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,  -3.0159372777109039e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 2.8348126021733977e-01, 9.1031840499614625e-01},
                {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 9.1417864987446651e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.5713389257830902e-01, 1.9165797370723034e-01}
        };
        //shifting
        for (int j = 0; j < xs.size(); j++) {
            x[0][j] = xs.get(j) - shiftedMatrix[j];
        }
        //rotating
        double shitedRotaatedX[][] = multiplyByMatrix(x, rotatedMatrix);


        for(int i =1; i <=D; i++){
            double xi = shitedRotaatedX[0][i-1];
            double inerSum_1 = 0.0;
            for(int k=0; k<=kMax; k++){
                inerSum_1 += Math.pow(a, k) * Math.cos(2 * Math.PI * Math.pow(b,k)*(xi+0.5));
            }
            sum += inerSum_1;
        }
        double inerSum_2 =0.0;
        for(int k=0; k<=kMax; k++){
            inerSum_2 += Math.pow(a, k) * Math.cos(Math.PI * Math.pow(b,k));
        }
        sum = sum - D * inerSum_2 +1;
        return sum +1;
    }

    public double CEC07(ArrayList<Double> xs){
        int D = xs.size();
        double zi = 0.0;
        double g =0.0;
        double x[][] = new double[1][D];
        double shiftedMatrix[] = {1.5519604466631876e+00,   3.7992270681072000e+00,   1.3609333677966774e+01,  -6.7928874412518397e+01,   7.9407748803220557e+01,   4.6034135728159043e+01,  -6.4280816830825444e+01,  -4.7688475683186425e+01,  -6.0210807314240753e+01,   3.6961469555721379e+01};
        double rotatedMatrix[][] = {
                {-3.4378315941460673e-02, -7.3911155710735865e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.0523399031716010e+00,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {1.1485242405257232e+00,  9.9172138327339543e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -9.8221173301823295e-01,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  8.6405702281889096e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -4.9170053952174497e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.0787048137178114e-01},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -3.7383444140566896e-01,  3.9203594526066760e-01,  0.0000000000000000e+00, -7.8796970160513635e-01,  0.0000000000000000e+00,  2.9267623305420509e-01,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  1.3610747065473336e-01, -6.0804319089793657e-01,  0.0000000000000000e+00, -5.6714757040795583e-01,  0.0000000000000000e+00, -5.3861105430076284e-01,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  8.5194811610211113e-02,  0.0000000000000000e+00,  0.0000000000000000e+00, -6.8358420399848130e-02,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  9.9401658458757081e-01},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -8.5125893429894595e-01,  8.0910637486172887e-03,  0.0000000000000000e+00,  2.3333420826786389e-01,  0.0000000000000000e+00, -4.6994458047268511e-01,  0.0000000000000000e+00},
                {1.2482031580671493e+00, -4.3601697061224165e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.7452756850803894e-01,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  3.4217070831253948e-01,  6.9030850372398578e-01,  0.0000000000000000e+00, -5.4795346377890963e-02,  0.0000000000000000e+00, -6.3513057403542883e-01,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  4.9613234664961658e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  8.6807701605010890e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  1.7175238382100888e-02}
        };
        //shifting
        for (int j = 0; j < xs.size(); j++) {
            x[0][j] = xs.get(j) - shiftedMatrix[j];
        }
        //rotating
        double shitedRotaatedX[][] = multiplyByMatrix(x, rotatedMatrix);

        for(int i =1; i <=D; i++){
            double xi = shitedRotaatedX[0][i-1];
            zi = xi+ 420.9687462275036;
            if(Math.abs(zi) <= 500) g += zi * Math.sin(Math.pow(Math.abs(zi), 1/2));
            else if(zi > 500) g += (500 - (zi % 500)) * Math.sin(Math.sqrt(Math.abs(500- (zi % 500))))- (Math.pow(zi - 500, 2)/(10000 * D));
            else if (zi < -500) g += ((Math.abs(zi) % 500) - 500)* Math.sin(Math.sqrt(Math.abs((zi % 500) - 500)))- (Math.pow(zi - 500, 2)/ (10000 * D));
        }

        return  (418.9829 * D - g) +1;
    }

    public double CEC08(ArrayList<Double> xs){
        int D = xs.size();
        double g = 0.0;

        double x[][] = new double[1][D];
        double shiftedMatrix[] = {7.5809536201790706e+01,   5.0874943496135501e+01,   1.5175339549395872e+01,   1.1931806696547099e+01,   5.7875148867198789e+01,   6.7627011010249618e+01,  -3.2825950734701912e+01,  -2.5753998135101980e+01,  -4.7446656658987820e+01,   4.0415323917015940e+00 };
        double rotatedMatrix[][] = {
                {3.2765524541169905e-01,   0.0000000000000000e+00,   9.4933157553264147e-01,   0.0000000000000000e+00,   0.0000000000000000e+00,   0.0000000000000000e+00,   0.0000000000000000e+00,  0.0000000000000000e+00,   0.0000000000000000e+00,  -6.2503273493416478e-01},
                {0.0000000000000000e+00, -3.9416689281102152e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  6.4866583791924892e-01,  0.0000000000000000e+00, -6.1484056671741649e-01,  2.1409383187478670e-01,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {8.9664357708578779e-01,  0.0000000000000000e+00,  5.3753076622704354e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  7.2208870049867158e-01},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.7334340479710361e-01,  0.0000000000000000e+00,  1.7833024809453213e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -9.6858163653240936e-01,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  5.5409637934619971e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -3.0250402339187688e-01,  0.0000000000000000e+00, -4.5590738714655138e-01,  6.2738901215463105e-01,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -9.6746517640794338e-01,  0.0000000000000000e+00,  1.5319756775156032e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  2.0134954102999447e-01,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  6.2594545303749249e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  6.9750321067293941e-01,  0.0000000000000000e+00,  3.4700267475621172e-01,  3.5646944254231774e-02,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  3.8184021897755271e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -3.4831274010388026e-02,  0.0000000000000000e+00, -5.4194896030516160e-01, -7.4784768097931698e-01,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.8429106449108246e-01,  0.0000000000000000e+00, -9.7197161884987071e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.4597251693095872e-01,  0.0000000000000000e+00},
                {1.7153665513641250e+00,  0.0000000000000000e+00, -7.8944707263781244e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -3.0034616211515719e-01}
        };
        //shifting
        for (int j = 0; j < xs.size(); j++) {
            x[0][j] = xs.get(j) - shiftedMatrix[j];
        }
        //rotating
        double shitedRotaatedX[][] = multiplyByMatrix(x, rotatedMatrix);


        for(int i=0, j =1; i<D-1 && j <D; i++, j++){
            double xi = shitedRotaatedX[0][i] ;
            double yi = shitedRotaatedX[0][j] ;
            g += 0.5 +( Math.pow(Math.sin(Math.sqrt(Math.pow(xi,2)+Math.pow(yi,2))), 2))- 0.5 / Math.pow((1+ 0.001 *(Math.pow(xi,2)+Math.pow(yi,2))), 2);
        }
        double xi = shitedRotaatedX[0][D-1] ;
        double yi = shitedRotaatedX[0][0] ;
        g += 0.5 +( Math.pow(Math.sin(Math.sqrt(Math.pow(xi,2)+Math.pow(yi,2))), 2))-0.5 / Math.pow((1+ 0.001 *(Math.pow(xi,2)+Math.pow(yi,2))), 2);
        return  g+1;
    }

    public double CEC09(ArrayList<Double> xs) {
        int D = xs.size();
        double sum1=0.0;
        double sum2=0.0;
        double sum3=0.0;
        double x[][] = new double[1][D];
        double shiftedMatrix[] = {-6.0107960952496171e+00,  -6.3449972860258995e+01,  -3.6938623728667750e+00,  -2.7449007717635965e+00,  -5.3547271030744199e+01,   3.1015786282259867e+01,   2.3200459416583499e+00,  -4.6987858548289097e+01,   3.5061378905112562e+01,  -3.4047417731046465e+00};
        double rotatedMatrix[][] = {
                {-7.6923624057192400e-02,  0.0000000000000000e+00,  7.2809258658661558e-02,  6.1371429917067155e-01,  0.0000000000000000e+00,  7.8239141541106805e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00, -1.1499983823069659e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  1.5729072158274271e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.3309066870600375e+00,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {-1.6730831752378217e-02,  0.0000000000000000e+00,  4.9480374519689890e-01,  6.5982384537901573e-01,  0.0000000000000000e+00, -5.6526261691115431e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {9.1421044415115027e-01,  0.0000000000000000e+00, -3.3249140365486585e-01,  2.2489758522716782e-01,  0.0000000000000000e+00, -5.5586027556918202e-02,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00, -1.2704704488967578e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -7.6341623484218024e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  4.9980922801223232e-01,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {-3.9751993551989051e-01,  0.0000000000000000e+00, -7.9957334378299227e-01,  3.7068629354440513e-01,  0.0000000000000000e+00, -2.5544478964007222e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  9.0184624725396623e-02,  0.0000000000000000e+00, -3.4243496198122719e-01, -9.3520320266563195e-01},
                {0.0000000000000000e+00,  9.7696981452382070e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -6.8376531090322690e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -4.7094671586086240e-01,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -6.1661688294908079e-01,  0.0000000000000000e+00, -7.5660067900409822e-01,  2.1757534831110126e-01},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -7.8208078427058847e-01,  0.0000000000000000e+00,  5.5704013261474550e-01, -2.7938492717261593e-01}
        };
        //shifting
        for (int j = 0; j < xs.size(); j++) {
            x[0][j] = xs.get(j) - shiftedMatrix[j];
        }
        //rotating
        double shitedRotaatedX[][] = multiplyByMatrix(x, rotatedMatrix);

        for(int i= 0; i <D; i++){
            double xi = shitedRotaatedX[0][i];
            sum1 +=Math.pow(xi, 2)- D;
            sum2 += Math.pow(xi, 2);
            sum3 += xi;
        }

        return  (Math.pow(Math.abs(sum1), 1/4)+(0.5 * sum2+sum3)/D + 0.5) +1;
    }

    public double CEC10(ArrayList<Double> xs){
        int D = xs.size();
        double sum1=0.0;
        double sum2=0.0;
        double x[][] = new double[1][D];
        double shiftedMatrix[] = {6.1441309549566370e-001,  1.8049534213689469e+001,  5.1107558757100151e+001,  5.1022671188681272e+000, -4.7667984552250942e+001, -7.3770454911164904e+000 ,-1.1534252828772665e+001 , 7.4568439937919834e+001,  1.9208808661355789e+001,  3.1262392306880571e+001};
        double rotatedMatrix[][] = {
                {-3.6144665808053256e-02,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.0275628429515489e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -9.9404965126067890e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  5.7732032557209267e-01,  6.3720045355332378e-01,  0.0000000000000000e+00, -1.2849837048835189e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  4.9413054191641514e-01,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00, -3.9586712917636224e-01, -2.8968310271397846e-01,  0.0000000000000000e+00, -5.1670886261404880e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  7.0170140895951061e-01,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {8.7364129527216294e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  4.7972086095511046e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -8.1355901812131204e-02,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  7.1410594914990388e-01, -6.7354832197676640e-01,  0.0000000000000000e+00, -1.9014210237853424e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.5209610582332084e-02,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -2.7911553415074425e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  1.6909147808954352e+00, -5.7964733149576031e-01},
                {-4.8522618471059614e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  8.7138340677473369e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -7.2432783108600157e-02,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00, -6.5689502748026013e-03,  2.3746987168001482e-01,  0.0000000000000000e+00, -8.2483095297218667e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -5.1304854346890183e-01,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.3268067876803948e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  7.1252154299165178e-04, -8.8929149979214617e-01},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -6.6352395365628714e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  4.6075742954200716e-01,  6.9770224192342367e-01}
        };
        //shifting
        for (int j = 0; j < xs.size(); j++) {
            x[0][j] = xs.get(j) - shiftedMatrix[j];
        }
        //rotating
        double shitedRotaatedX[][] = multiplyByMatrix(x, rotatedMatrix);

        for(int i =0 ; i<D; i++){
            double xi = shitedRotaatedX[0][i];
            sum1 += Math.pow(xi, 2);
            sum2 += Math.cos(2 * Math.PI*xi);
        }

        return (-20* Math.exp((0.2*Math.sqrt(1/D * sum1))) - Math.exp(1/D * sum2)+20+Math.E)+1;
    }

    int getDimensions() {
        int dimensions = 2;
        switch (functionName) {
            case "function_FM":
                dimensions = 6;
                break;
            case "function_antenna":
                dimensions = 4;
                break;
            case "FX1":
                dimensions = 10;
                break;
            case "FX2":
                dimensions = 10;
                break;
            case "FX3":
                dimensions = 10;
                break;
            case "FX4":
                dimensions = 10;
                break;
            case "FX5":
                dimensions = 10;
                break;
            case "FX6":
                dimensions = 10;
                break;
            case "FX7":
                dimensions = 10;
                break;
            case "FX8":
                dimensions = 10;
                break;
            case "FX9":
                dimensions = 10;
                break;
            case "FX10":
                dimensions = 10;
                break;
            case "FX11":
                dimensions = 10;
                break;
            case "FX12":
                dimensions = 10;
                break;
            case "FX13":
                dimensions = 10;
                break;
            case "FX14":
                dimensions = 2;
                break;
            case "FX15":
                dimensions = 4;
                break;
            case "FX16":
                dimensions = 10;
                break;
            case "FX17":
                dimensions = 10;
                break;
            case "FX18":
                dimensions = 10;
                break;
            case "FX19":
                dimensions = 10;
                break;
            case "FX20":
                dimensions = 5;
                break;
            case "FX21":
                dimensions = 2;
                break;
            case "FX22":
                dimensions = 5;
                break;
            case "FX23":
                dimensions = 5;
                break;
            case "FX24":
                dimensions = 5;
                break;
            case "CEC01":
                dimensions = 9;
                break;
            case "CEC02":
                dimensions = 16;
                break;
            case "CEC03":
                dimensions = 18;
                break;
            case "CEC04":
            case "CEC05":
            case "CEC06":
            case "CEC07":
            case "CEC08":
            case "CEC09":
            case "CEC10":
                dimensions = 10;
                break;
        }
        return dimensions;
    }

    double[] getUpperBound() {
        double[] upperBound = {0, 0};
        switch (functionName) {
            case "function_FM":
                upperBound[0] = 6.35;
                break;
            case "function_antenna":
                upperBound[0] = 2.0;
                break;
            case "FX1":
                upperBound[0] = 100;
                break;
            case "FX2":
                upperBound[0] = 10;
                break;
            case "FX3":
                upperBound[0] = 100;
                break;
            case "FX4":
                upperBound[0] = 100;
                break;
            case "FX5":
                upperBound[0] = 30;
                break;
            case "FX6":
                upperBound[0] = 100;
                break;
            case "FX7":
                upperBound[0] = 1.28;
                break;
            case "FX8":
                upperBound[0] = 500;
                break;
            case "FX9":
                upperBound[0] = 5.12;
                break;
            case "FX10":
                upperBound[0] = 32;
                break;
            case "FX11":
                upperBound[0] = 600;
                break;
            case "FX12":
                upperBound[0] = 50;
                break;
            case "FX13":
                upperBound[0] = 50;
                break;
            case "FX14":
                upperBound[0] = 65.536;
                break;
            case "FX15":
                upperBound[0] = 5;
                break;
            case "FX16":
                upperBound[0] = 5;
                break;
            case "FX17":
                upperBound[0] = 5;
                break;
            case "FX18":
                upperBound[0] = 5;
                break;
            case "FX19":
                upperBound[0] = 5;
                break;
            case "FX20":
                upperBound[0] = 1;
                break;
            case "FX21":
                upperBound[0] = 1;
                break;
            case "FX22":
                upperBound[0] = 1;
                break;
            case "FX23":
                upperBound[0] = 1;
                break;
            case "FX24":
                upperBound[0] = 1;
                break;
            case "CEC01":
                upperBound[0] = 8192;
                break;
            case "CEC02":
                upperBound[0] = 16384;
                break;
            case "CEC03":
                upperBound[0] = 4;
                break;
            case "CEC04":
            case "CEC05":
            case "CEC06":
            case "CEC07":
            case "CEC08":
            case "CEC09":
            case "CEC10":
                upperBound[0] = 100;
                break;
        }
        return upperBound;
    }

    double[] getLowerBound() {
        double[] Lower = {0, 0};
        switch (functionName) {

            case "function_FM":
                Lower[0] = -6.4;
                break;
            case "function_antenna":
                Lower[0] = 0;
                break;
            case "FX1":
                Lower[0] = -100;
                break;
            case "FX2":
                Lower[0] = -10;
                break;
            case "FX3":
                Lower[0] = -100;
                break;
            case "FX4":
                Lower[0] = -100;
                break;
            case "FX5":
                Lower[0] = -30;
                break;
            case "FX6":
                Lower[0] = -100;
                break;
            case "FX7":
                Lower[0] = -1.28;
                break;
            case "FX8":
                Lower[0] = -500;
                break;
            case "FX9":
                Lower[0] = -5.12;
                break;
            case "FX10":
                Lower[0] = -32;
                break;
            case "FX11":
                Lower[0] = -600;
                break;
            case "FX12":
                Lower[0] = -50;
                break;
            case "FX13":
                Lower[0] = -50;
                break;
            case "FX14":
                Lower[0] = -65.536;
                break;
            case "FX15":
                Lower[0] = -5;
                break;
            case "FX16":
                Lower[0] = -5;
                break;
            case "FX17":
                Lower[0] = -5;
                break;
            case "FX18":
                Lower[0] = -5;
                break;
            case "FX19":
                Lower[0] = -5;
                break;
            case "FX20":
                Lower[0] = 0;
                break;
            case "FX21":
                Lower[0] = 0;
                break;
            case "FX22":
                Lower[0] = 0;
                break;
            case "FX23":
                Lower[0] = 0;
                break;
            case "FX24":
                Lower[0] = 0;
                break;
            case "CEC01":
                Lower[0] = -8192;
                break;
            case "CEC02":
                Lower[0] = -16384;
                break;
            case "CEC03":
                Lower[0] = -4;
                break;
            case "CEC04":
            case "CEC05":
            case "CEC06":
            case "CEC07":
            case "CEC08":
            case "CEC09":
            case "CEC10":
                Lower[0] = -100;
                break;
        }
        return Lower;
    }

    public double[][] multiplyByMatrix(double[][] m1, double[][] m2) {
        int m1ColLength = m1[0].length; // m1 columns length
        int m2RowLength = m2.length;    // m2 rows length
        if (m1ColLength != m2RowLength) return null; // matrix multiplication is not possible
        int mRRowLength = m1.length;    // m result rows length
        int mRColLength = m2[0].length; // m result columns length
        double[][] mResult = new double[mRRowLength][mRColLength];
        for (int i = 0; i < mRRowLength; i++) {         // rows from m1
            for (int j = 0; j < mRColLength; j++) {     // columns from m2
                for (int k = 0; k < m1ColLength; k++) { // columns from m1
                    mResult[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
        return mResult;
    }


}