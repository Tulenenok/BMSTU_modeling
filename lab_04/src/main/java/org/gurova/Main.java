package org.gurova;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

public class Main {
    public static void main(String[] args) throws Exception {
        Model mod = new Model(3, 1.4, 0.35, 0.5, 1, 300, 5.668e-12, 100.0, 0.05, 2e-3);
        int n = ((int)((mod.zN - mod.z0) / mod.h / mod.R)) + 1;
        double[] z = new double[n];
        double[] y = new double[n];
        double[] x = new double[n];
        double[] t = new double[n];

        List<double[]> ys = new ArrayList<>();

        for (int i = 0; i < n; i++)
        {
            z[i] = mod.z0 + i * mod.h / mod.R;
            x[i] = mod.z0 + i * mod.h * mod.R;
            y[i] = mod.T0;
        }

        for (double i = 0; i < 200; i += mod.tau)
        {
            double[] tmp = new double[n];
            System.arraycopy(y, 0, tmp, 0, y.length);           // y.CopyTo(tmp, 0);
            ys.add(tmp);
            y = mod.NextTime(z, y, i);
        }

        Color[] colors = {Color.BLUE};

        /* Граф 1 про Fo */
        List<Double> lx = new ArrayList<>();
        List<Double> ly = new ArrayList<>();

        for(double i = 0; i < 1000; i++) {
            lx.add(i);
            ly.add(mod.F0(i));
        }

        new Graph("ЛР №4", "График изменения Fo", lx, ly, "F0", "Время t", "F0", colors);
        /**/


        /* Граф 2 про температуру на внутренней стенке трубки */
        lx = new ArrayList<>();
        ly = new ArrayList<>();

        for(int i = 0; i < ys.size(); i++) {
            lx.add(i * mod.tau);
            ly.add(ys.get(i)[0]);
        }

        new Graph("ЛР №4", "График изменения температуры на внутренней стенке трубы", lx, ly, "",
                "Время t", "Температура", colors);
        /**/


        /* Граф 3 про много графиков с температурой */
        Graph graph = new Graph("ЛР №4", "График изменения температуры по удалению от стенки трубы");
        List<String> names = new ArrayList<>();
        for(int i = 0; i < 15; i++) {
            names.add("T " + i);
        }

        graph.addData(ys, x, names, n);
        graph.fillChart("Ось х", "Ось y", null);
        graph.show();
        /**/
    }
}