package com.github.jasonmerecki.optpricing.black;

import static  com.github.jasonmerecki.optpricing.util.NormalDistribution.StandardNormal.cdf;

// references:
// http://www.macroption.com/black-scholes-formula/
// https://introcs.cs.princeton.edu/java/22library/BlackScholes.java.html

public class BlackLike {

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
		
		double re = Math.exp(-r*t), qe = Math.exp(-q*t);
		double d1 = d1(s, k, t, v, r, q);
		double d2 = d2(d1, v, t);
		double nd1, nd2;
		
		d1 = sign * d1; 
		d2 = sign * d2;
		nd1 = cdf(d1);
		nd2 = cdf(d2);
		
		double bsprice = sign * ((s * qe * nd1) - (k * re * nd2));
		return bsprice;
	}
	
	private static double d1(double s, double k, double t, double v, double r, double q) {
		double d1 = Math.log(s/k) +  (t * (r - q + ((v * v) * 0.5d) ));
		d1 = d1 / (v * (Math.sqrt(t)));
		return d1;
	}
	
	private static double d2(double d1, double v, double t) {
		double vt = (v * (Math.sqrt(t)));
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
		// NEED TO CHECK: in the money American options are priced at intrinsic value
		if (k < s) {
			// return Math.abs(s - k);
		}
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

}
