package com.github.jasonmerecki.optpricing.black;

import static  com.github.jasonmerecki.optpricing.util.NormalDistribution.StandardNormal.*;

// references:
// http://www.macroption.com/black-scholes-formula/
// https://introcs.cs.princeton.edu/java/22library/BlackScholes.java.html

public class BlackLike {
	
	public static final double IV_PRECISION = 0.00001d;
	private static final double ONE_DAY = 1d / 365d;
	private static final double SQ_TWOPI = (Math.sqrt(2 * Math.PI));

	// type == C(all) or P(ut)
	// s = stock price (current), k = strike price
	// t = expiry time (annualized where 1 = one year), v = volatility, r = risk-free rate
	// q = dividend yield
	
	public static double priceBlackScholes(String type, double s, double k, double t, double v, double r, double q) {
		double sign = 0d;
		if ("C".equals(type)) {
			if (t <= 0l) {
				return Math.abs(s - k);
			}
			sign = 1;
		} else if ("P".equals(type)) {
			if (t <= 0l) {
    				return Math.abs(k - s);
			}
			sign = -1;
		} else {
			return 0.0d;
		}
		
		double dr = Math.exp(-r*t), drq = Math.exp(-q*t);
		double vt = (v * (Math.sqrt(t)));
		double d1 = d1(s, k, t, v, r, q, vt);
		double d2 = d2(d1, vt);
		double nd1, nd2;
		
		d1 = sign * d1; 
		d2 = sign * d2;
		nd1 = cdf(d1);
		nd2 = cdf(d2);
		
		double bsprice = sign * ((s * drq * nd1) - (k * dr * nd2));
		return bsprice;
	}
	
	private static double d1(double s, double k, double t, double v, double r, double q, double vt) {
		double d1 = Math.log(s/k) +  (t * (r - q + ((v * v) * 0.5d) ));
		d1 = d1 / vt;
		return d1;
	}
	
	private static double d2(double d1, double vt) {
		double d2 = d1 - vt;
		return d2;
	}
	
	public static double priceBjerkStens(String type, double s, double k, double t, double v, double r, double q) {
	    	double b = r - q; // b is the cost of carry
	    	double price = 0.0d;
	    if ("P".equals(type)) {
	    		if (t <= 0l) {
	    			return Math.abs(k - s);
	    		}
	    		if (-b > q || q == 0) {
	    			return priceBlackScholes(type, s, k, t, v, r, q);
	    		}
	    		price = bjerkStensInternal(type, k, s, t, v, -b, r - b, q);
	    } else if ("C".equals(type)) {
	    		if (t <= 0l) {
				return Math.abs(s - k);
			}
	    		if (b > r || q == 0) {
	    			return priceBlackScholes(type, s, k, t, v, r, q);
	    		}
	    		price = bjerkStensInternal(type, s, k, t, v, b, r, q);
	    }
		return price;
	}
	
	private static double bjerkStensInternal(String type, double s, double k, double t, double v, double b, double r, double q) {
		double v2 = v * v;
		double beta = (0.5 - (b / v2)) + Math.pow( ( Math.pow(( (b / v2) - 0.5), 2) + (2 * r / v2) ), 0.5);

	    	double betainfinity = ( beta / (beta -1) ) * k;	
	    	double betazero = Math.max(k, (r / q) * k);
	    	
	    	double h = -((b * t) + (2 * v * Math.pow(t, 0.5))) * ( (k * k) / ((betainfinity - betazero) * betazero)); 
	    	double X = betazero + ((betainfinity - betazero) * (1 - Math.exp(h)));
	    	
	    	if (X < s) {
	    		// price equals intrinsic value
	    		return Math.abs(s - k);
	    	} else {
	    		double alpha = (X -k) * Math.pow(X, -beta);
	    		double tmp1 = alpha * Math.pow( s, beta);
	    		double tmp2 = alpha * bjerkStensPhi(s, t, beta, X, X, v, r, b);
	    		double tmp3 = bjerkStensPhi(s, t, 1, X, X, v, r, b);
	    		double tmp4 = bjerkStensPhi(s, t, 1, k, X, v, r, b);
	    		double tmp5 = k * bjerkStensPhi(s, t, 0, X, X, v, r, b);
	    		double tmp6 = k * bjerkStensPhi(s, t, 0, k, X, v, r, b);
	    		double bjprice = tmp1 - tmp2 + tmp3 - tmp4 - tmp5 + tmp6;
	    		return bjprice;
	    	}
	}
	
	private static double bjerkStensPhi(double s, double t, double gamma, double h, double X, double v, double r, double b) {
		double v2 = v * v;
	    	double K = ((2 * b) / v2) + (2 * gamma - 1);
	    	double lambda = -r + (gamma * b) + ((0.5 * gamma * (gamma - 1)) * v2);
	    	double tmp1 = ((Math.log(s / h)) + (b + (gamma - 0.5) * v2) * t) / (v * Math.pow(t, 0.5));
	    	double tmp2 = (Math.log((X * X) / (s * h)) + (b + (gamma - 0.5) * v2) * t) / (v * Math.pow(t, 0.5)); 
	    	tmp2 = tmp1 + 2 * Math.log(X / s) / (v * Math.pow(t, 0.5));
	    	
	    	double psphi = Math.exp(lambda * t) * Math.pow(s, gamma) * (cdf(-tmp1) - (Math.pow(X / s, K)) * cdf(-tmp2));
		return psphi;
	}
	
	private static double d1pdf(double s, double k, double v,
			double t, double r, double q) {
	    double vt = (v * (Math.sqrt(t)));
	    double d1 = d1(s, k, t, v, r, q, vt);
	    double etod1sqhalf = Math.exp(-(d1*d1) * 0.5);
	    etod1sqhalf = etod1sqhalf/SQ_TWOPI;
	    return etod1sqhalf;
	}
	
	// Greeks (generally follows the macroption.com spreadsheet formula)
	public static double bsVega (String type, double s, double k, double v,
			double t, double r, double q) {
		double d1pdf = d1pdf(s, k, v, t, r, q);
		double drq = Math.exp(-q*t);
		double sqt = Math.sqrt(t);
		double vega = (d1pdf) * drq * s * sqt * 0.01;
		return vega;
	}
	
	public static double bjVega (String type, double s, double k, double v,
			double t, double r, double q) {
		double m = 0.01;
		double val1 = priceBjerkStens(type, s, k, t, v, r, q);
		double val2 = priceBjerkStens(type, s, k, t, v + m, r, q);
		double vega = (val2 - val1) / m / 100d;
		return vega;	    	    
	}
	
	public static double bsDelta (String type, double s, double k, double v,
			double t, double r, double q) {
	    double drq = Math.exp(-q*t);
	    double zo = ("P".equals(type)) ? -1d : 0d;
	    double vt = (v * (Math.sqrt(t)));
	    double d1 = d1(s, k, t, v, r, q, vt);
	    double cdfd1 = cdf(d1);
	    double delta = drq * (cdfd1 - zo);
	    return delta;
	}
	
	public static double bjDelta (String type, double s, double k, double v,
			double t, double r, double q) {
		double umove = 1.01;
		double dmove = 1 / umove;
		double uval = priceBjerkStens(type, s * umove, k, t, v, r, q);
		double dval = priceBjerkStens(type, s * dmove, k, t, v, r, q);
		double delta = (uval - dval) / (s * (umove - dmove));
		return delta;
	}
	
	public static double bsGamma (String type, double s, double k, double v,
			double t, double r, double q) {
		double drq = Math.exp(-q*t);
		double drd = (s * v * Math.pow(t, 0.5));
		double d1pdf = d1pdf(s, k, v, t, r, q);
		double gamma = (drq / drd)  * d1pdf;
	    return gamma;
	}
	
	public static double bsTheta (String type, double s, double k, double v,
			double t, double r, double q) {
		double sign = ("P".equals(type)) ? -1d : 1d;
		double drq = Math.exp(-q*t);
		double dr = Math.exp(-r*t);
		double d1pdf = d1pdf(s, k, v, t, r, q);
		double twosqt = 2 * Math.sqrt(t);
		double p1 = -1 * ((s * v * drq)/twosqt) * d1pdf;

		double vt = (v * (Math.sqrt(t)));
		double d1 = d1(s, k, t, v, r, q, vt);
		double d2 = d2(d1, vt);
		double nd1, nd2;
		
		d1 = sign * d1; 
		d2 = sign * d2;
		nd1 = cdf(d1);
		nd2 = cdf(d2);
		
		double p2 = -sign * r * k * dr * nd2;
		double p3 = sign * q * s * drq * nd1;
		double theta = (p1 + p2 + p3) / 365;
		return theta;
	}
	
	public static double bjTheta (String type, double s, double k, double v,
			double t, double r, double q) {
		double y = t - ONE_DAY;
		y = (y < 0) ? 0 : y;
		double val = priceBjerkStens(type, s, k, t, v, r, q);
		double valt =  priceBjerkStens(type, s, k, y, v, r, q);
		return valt - val;
	}
	
	public static double bsRho (String type, double s, double k, double v,
			double t, double r, double q) {
		double sign = ("P".equals(type)) ? -1d : 1d;
		double dr = Math.exp(-r*t);
		double p1 = sign * (k * t * dr) / 100;
		
		double vt = (v * (Math.sqrt(t)));
		double d1 = d1(s, k, t, v, r, q, vt);
		double d2 = sign * d2(d1, vt);
		double nd2 = cdf(d2);
		double rho = p1 * nd2;
		return rho;
	}

	// Implied vol
	public static double bsImpliedVol(String type, double p, double s, 
			double k, double r, double t, double v, double q) {
		v = v == 0d ? 0.5 : v;
		double errlimit = IV_PRECISION;
		double maxloops = 100;
		double dv = errlimit + 1;
		double n = 0;
		while (Math.abs(dv) > errlimit && n < maxloops) {
			double difval = priceBlackScholes(type, s, k, t, v, r, q) - p;
			double v1 = bsVega(type, s, k, v, t, r, q) / 0.01;
			dv = difval / v1;
			v = v - dv;
			n++;
		}
		return n < maxloops ? v : Double.NaN;
	}
	
	public static double bjImpliedVol(String type, double p, double s, 
			double k, double r, double t, double v, double q) {
		v = v == 0d ? 0.5 : v;
		double errlimit = IV_PRECISION;
		double maxloops = 100;
		double dv = errlimit + 1;
		double n = 0;
		while (Math.abs(dv) > errlimit && n < maxloops) {
			double difval = priceBjerkStens(type, s, k, t, v, r, q) - p;
			double v1 = bjVega(type, s, k, v, t, r, q) / 0.01;
			dv = difval / v1;
			v = v - dv;
			n++;
		}
		return n < maxloops ? v : Double.NaN;
	}

}
