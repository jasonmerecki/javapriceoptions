package com.github.jasonmerecki.optpricing.black;

import static com.github.jasonmerecki.optpricing.util.NormalDistribution.*;

// references:
// http://www.macroption.com/black-scholes-formula/

public class BlackLike {

	// type == C(all) or P(ut)
	// s = stock price (current), k = strike price
	// t = expiry time, v = volatility, r = risk-free rate
	// q = dividend yield
	
	public static double BlackScholes(String type, double s, double k, double t, double v, double r, double q) {
		
		double re = Math.exp(-r*t), qe = Math.exp(-q*t);
		double d1 = d1(s, k, t, v, r, q);
		double d2 = d2(d1, v, t);
		double nd1, nd2;
		double sign = 0d;
		if ("C".equals(type)) {
			sign = 1;
		} else if ("P".equals(type)) {
			sign = -1;
		} else {
			return 0.0d;
		}
		d1 = sign * d1; 
		d2 = sign * d2;
		nd1 = StandardNormal.cdf(d1);
		nd2 = StandardNormal.cdf(d2);
		
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

}
