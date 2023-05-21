package org.gurova;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.awt.*;
import java.util.List;

public class Graph {
    private JFreeChart chart;
    private String bigTitle;
    private String smallTitle;
    private XYSeriesCollection dataset;

    public Graph(String bigTitle, String smallTitle) {
        this.bigTitle = bigTitle;
        this.smallTitle = smallTitle;
        this.dataset = new XYSeriesCollection();
    }

    public Graph(String bigTitle, String smallTitle, List<Double> x, List<Double> y, String name,
                 String xLabel, String yLabel, Color[] colors) throws Exception {

        this.bigTitle = bigTitle;
        this.smallTitle = smallTitle;
        this.dataset = new XYSeriesCollection();

        this.addData(x, y, name);
        fillChart(xLabel, yLabel, colors);
        show();
    }

    public void addData(List<Double> x, List<Double> y, String name) throws Exception {
        if (x.size() != y.size()) {
            throw new Exception("Размеры массивов х и y не совпадают");
        }

        XYSeries series = new XYSeries(name);

        for (int i = 0; i < x.size(); i++) {
            series.add(x.get(i), y.get(i));
        }

        dataset.addSeries(series);
    }

    public void fillChart(String xLabel, String yLabel, Color[] colors) {
        chart = ChartFactory.createXYLineChart(
                smallTitle,
                xLabel,
                yLabel,
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        if (colors != null) {
            XYPlot plot = chart.getXYPlot();

            for(int i = 0; i < colors.length; i++) {
                plot.getRenderer().setSeriesPaint(i, colors[i]);
            }
        }
    }

    public void show() {
        ChartFrame frame = new ChartFrame(bigTitle, chart);
        frame.pack();
        frame.setVisible(true);
    }

    public void addData(List<double[]> ys, double[] x,  List<String> names, int n) {
        // Закостыленный метод чисто для этой лабы !!!

        for (int j = 0; j < 15; j++) {
            XYSeries series = new XYSeries(names.get(j));

            for (int i = 0; i < n; i++) {
                series.add(x[i], ys.get((ys.size() - 1) / 15 * j)[i]);
            }

            dataset.addSeries(series);
        }

    }
}
