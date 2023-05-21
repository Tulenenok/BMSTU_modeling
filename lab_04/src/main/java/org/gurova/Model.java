package org.gurova;

import java.util.Arrays;

public class Model {
    public double Np;
    public double r0;
    public double R;
    public double z0;
    public double zN;
    public double T0;
    public double sigma;
    public double alpha;
    public double h;
    public double a2 = 2.049;
    public double b2 = 0.563e-3;
    public double c2 = 0.528e5;
    public double m2 = 1;
    public double tau = 4e-6;
    public double Fmax = 5000;
    public double tmax = 150e-6;

    public Model()
    {
        z0 = r0 / R;
    }

    public Model(double tau, double np, double r0, double r, double zN, double t0, double sigma, double f0, double alpha, double h)
    {
        this.tau = tau;
        Np = np;
        this.r0 = r0;
        R = r;
        this.zN = zN;
        T0 = t0;
        this.sigma = sigma;
        this.alpha = alpha;
        this.h = h;
        z0 = r0 / R;
    }

    public static double Labmda(double T)
    {
        //return 2.50e-2;
        double[] lambda = new double[] { 1.36e-2, 1.63e-2, 1.81e-2, 1.98e-2, 2.50e-2, 2.74e-2 };
        double[] t = new double[]      {     300,     500,     800,    1100,     2000,   2400 };
        return Interpolation(T, t, lambda);
    }

    private static double K(double T)
    {
        double[] k = new double[] { 2.0e-2, 5.0e-2, 7.8e-2, 1.0e-1, 1.3e-1, 2.0e-1 };
        double[] t = new double[] { 293, 1278, 1528, 1677, 2000, 2400 };
        return Interpolation(T, t, k);
    }

    public double F0(double t)
    {
        if (t < 100) {
            return 100;
        }
        else {
            return 0;
        }

//        return Fmax / tmax * t * Math.exp(-(t / tmax - 1));
    }

    private double P(double T)
    {
        return 4 * Np * Np * sigma * K(T) * Math.pow(T, 3);
    }

    private double F(double T)
    {
        return 4 * Np * Np * sigma * K(T) * Math.pow(T0, 4);
    }

    private double V(double zn_1, double zn, double zn1)
    {
        return (Math.pow((zn1 + zn) / 2, 2) - Math.pow((zn + zn_1) / 2, 2)) / 2;
    }

    private double Ct(double yn)
    {
        return a2 + b2 * Math.pow(yn, m2) - c2 / (yn * yn);
    }

    private double A(double yn_1, double yn, double zn_1, double zn)
    {
        double ln_12 = (Labmda(yn_1) + Labmda(yn)) / 2;
        double zn_12 = (zn_1 + zn) / 2;
        return zn_12 * ln_12 / (R * R * h) * tau;
    }

    private double B(double yn_1, double yn, double yn1, double zn_1, double zn, double zn1)
    {
        return A(yn_1, yn, zn_1, zn) + C(yn, yn1, zn, zn1) + P(yn) * V(zn_1, zn, zn1) * tau + zn * Ct(yn) * h;
    }

    private double C(double yn, double yn1, double zn, double zn1)
    {
        double ln12 = (Labmda(yn1) + Labmda(yn)) / 2;
        double zn12 = (zn1 + zn) / 2;
        return zn12 * ln12 / (R * R * h) * tau;
    }

    private double D(double yn, double zn_1, double zn, double zn1, double ym)
    {
        return F(yn) * V(zn_1, zn, zn1) * tau + zn * Ct(yn) * h * ym;
    }

    static double Interpolation(double xval, double[] xs, double[] ys)
    {
        int i = 0;
        while (i < xs.length && xval < xs[i]) {
            i++;
        }

        if (i >= xs.length - 1) {
            i = xs.length - 2;
        }

        return ys[i] + (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]) * (xval - xs[i]);
    }

    public double[] RunMethod(double[] a, double[] b, double[] c, double[] d, double m0, double k0, double mn, double kn, double p0, double pn)
    {
        double[] y = new double[b.length];
        double[] xi = new double[b.length - 1];
        double[] eta = new double[b.length - 1];

        xi[0] = -k0 / m0;
        eta[0] = p0 / m0;

        for (int i = 1; i < xi.length; i++)
        {
            double tmp = b[i] - a[i] * xi[i - 1];
            xi[i] = c[i] / tmp;
            eta[i] = (a[i] * eta[i - 1] + d[i]) / tmp;
        }

        y[y.length - 1] = (pn - mn * eta[eta.length - 1]) / (kn + mn * xi[xi.length - 1]);

        for (int i = y.length - 1; i > 0; i--) {
            y[i - 1] = xi[i - 1] * y[i] + eta[i - 1];
        }

        return y;
    }

    private double P0(double z0, double z1, double y0, double y1, double t, double ym0, double ym1)
    {
        double z12 = (z0 + z1) / 2;
        double f12 = (F(y0) + F(y1)) / 2;
        double c12 = (Ct(y0) + Ct(y1)) / 2;
        double ym12 = (ym0 + ym1) / 2;

        return (z0 * Ct(y0) * ym0 + z12 * c12 * ym12) * h / 2
                + z0 / R * F0(t) * tau + h / 4 * (F(y0) * z0 + f12 * z12) * tau;

        //return F0;
    }

    private double K0(double z0, double z1, double y0, double y1)
    {
        double z12 = (z0 + z1) / 2;
        double x12 = (Labmda(y0) + Labmda(y1)) / 2;
        double p12 = (P(y0) + P(y1)) / 2;
        double c12 = (Ct(y0) + Ct(y1)) / 2;

        return h / 4 * c12 * z12 - z12 * x12 * tau / (R * R * h) + h / 8 * p12 * z12 * tau;
        //return -x12 / (R * h);
    }

    private double M0(double z0, double z1, double y0, double y1)
    {
        double z12 = (z0 + z1) / 2;
        double x12 = (Labmda(y0) + Labmda(y1)) / 2;
        double p12 = (P(y0) + P(y1)) / 2;
        double c12 = (Ct(y0) + Ct(y1)) / 2;

        return h / 2 * (c12 * z12 / 2 + z0 * Ct(y0)) + z12 * x12 * tau / (R * R * h)
                + h * tau / 4 * (p12 * z12 / 2 + P(y0) * z0);

        // return x12 / (R * h);
    }

    private double PN(double zN_1, double zN, double yN_1, double yN, double ymN, double ymN_1)
    {
        double zN_12 = (zN_1 + zN) / 2;
        double fN_12 = (F(yN_1) + F(yN)) / 2;
        double cN_12 = (Ct(yN) + Ct(yN_1)) / 2;

        return (zN_12 * cN_12 + zN * Ct(yN)) * h / 2 + alpha * T0 * zN * tau / R
                + h / 4 * (F(yN) * zN + fN_12 * zN_12) * tau;

        // return alpha * T0;
    }

    private double KN(double zN_1, double zN, double yN_1, double yN)
    {
        double zN_12 = (zN + zN_1) / 2;
        double xN_12 = (Labmda(yN_1) + Labmda(yN)) / 2;
        double pN_12 = (P(yN) + P(yN_1)) / 2;
        double cN_12 = (Ct(yN) + Ct(yN_1)) / 2;

        return h / 2 * (zN_12 * cN_12 / 2 + zN * Ct(yN)) + zN_12 * xN_12 * tau / (R * R * h)
                + alpha * zN * tau / R + h * tau / 4 * (pN_12 * zN_12 / 2 + P(yN) * zN);

        // return -xN_12 / (R * h) + alpha;
    }

    private double MN(double zN_1, double zN, double yN_1, double yN)
    {
        double zN_12 = (zN + zN_1) / 2;
        double xN_12 = (Labmda(yN) + Labmda(yN_1)) / 2;
        double pN_12 = (P(yN) + P(yN_1)) / 2;
        double cN_12 = (Ct(yN) + Ct(yN_1)) / 2;

        return h / 4 * zN_12 * cN_12 - zN_12 * xN_12 * tau / (R * R * h) + h / 8 * pN_12 * zN_12 * tau;

        //return xN_12 / (R * h);
    }

    public double[] NextTime(double[] z, double[] ym, double t)
    {
        int n = z.length;
        double[] y = new double[n];

        System.arraycopy(ym, 0, y, 0, ym.length);           // ym.CopyTo(y, 0);

        double[] erry = new double[n];
        double maxerr;
        int its = 0;

        do
        {
            double[] a = new double[n], b = new double[n], c = new double[n], d = new double[n];
            double k0, m0, p0, kn, pn, mn;
            for (int i = 1; i < n - 1; i++)
            {
                a[i] = A(y[i - 1], y[i], z[i - 1], z[i]);
                b[i] = B(y[i - 1], y[i], y[i + 1], z[i - 1], z[i], z[i + 1]);
                c[i] = C(y[i], y[i + 1], z[i], z[i + 1]);
                d[i] = D(y[i], z[i - 1], z[i], z[i + 1], ym[i]);
            }
            k0 = K0(z[0], z[1], y[0], y[1]);
            m0 = M0(z[0], z[1], y[0], y[1]);
            p0 = P0(z[0], z[1], y[0], y[1], t, ym[0], ym[1]);
            kn = KN(z[n - 2], z[n - 1], y[n - 2], y[n - 1]);
            mn = MN(z[n - 2], z[n - 1], y[n - 2], y[n - 1]);
            pn = PN(z[n - 2], z[n - 1], y[n - 2], y[n - 1], ym[n - 1], ym[n - 2]);
            double[] newy = RunMethod(a, b, c, d, m0, k0, mn, kn, p0, pn);

            for (int i = 0; i < n; i++)
                erry[i] = Math.abs((newy[i] - y[i]) / newy[i]);
            y = newy;

            maxerr = Arrays.stream(erry).summaryStatistics().getMax();      // erry.Max();
            its++;
        }

        while (maxerr > 1e-3 &&  its < 50);
        System.out.println(its);
        return y;
    }
}
