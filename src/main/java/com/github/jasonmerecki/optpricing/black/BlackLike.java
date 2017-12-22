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
		double p1, nd1, nd2, p2;
		if ("C".equals(type)) {
			nd1 = StandardNormal.cdf(d1);
			nd2 = StandardNormal.cdf(d2);
			p1 = (s * qe) * nd1;
			p2 = (k * qe) * nd2;
		} else if ("P".equals(type)) {
			nd1 = StandardNormal.cdf(-d2);
			nd2 = StandardNormal.cdf(-d1);
			p1 = (k * re) * nd1;
			p2 = (s * qe) * nd2;
		} else {
			return 0.0d;
		}
		double bsprice = p1 - p2;
		return bsprice;
	}
	
	private static double d1(double s, double k, double t, double v, double r, double q) {
		double d1 = Math.log(s/k) +  (t * (r - q + ((v * v) * 0.5d) ));
		d1 = d1 / (v * (Math.sqrt(t)));
		return d1;
	}
	
	private static double d2(double d1, double v, double t) {
		double wtf = Math.sqrt(t);
		double vt = (v * (Math.sqrt(t)));
		double d2 = d1 - vt;
		return d2;
	}

}
