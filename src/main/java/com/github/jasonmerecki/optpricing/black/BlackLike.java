package com.github.jasonmerecki.optpricing.black;

import static  com.github.jasonmerecki.optpricing.util.NormalDistribution.StandardNormal.cdf;

// references:
// http://www.macroption.com/black-scholes-formula/

public class BlackLike {

	// type == C(all) or P(ut)
	// s = stock price (current), k = strike price
	// t = expiry time (annualized where 1 = one year), v = volatility, r = risk-free rate
	// q = dividend yield
	
	public static double BlackScholes(String type, double s, double k, double t, double v, double r, double q) {
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
	
	public static double BjerkStensPrice(String type, double s, double k, double t, double v, double r, double q) {
		
		/*
		 * Function BJValue(CallPutFlag, underlyerType, S, K, r, Dividend, sigma, t)
'Based on paper CLOSED FORM VALUATION OF AMERICAN OPTIONS
' http://www.nhh.no/for/dp/2002/0902.pdf
'CallPutFlag = p or c for put/call
'S = spot price
'K = strike
'r=risk free rate
'dividend = dividend yield
'sigma= volatility
'T=option maturity
*/

	    	double b = r - q; // b is the cost of carry
	    	double price = 0.0d;
	    // If (CallPutFlag = "p" Or CallPutFlag = "P") Then
	    if ("P".equals(type)) {
	        // 'use put-call transformation as give in eq (9)
	        // If t <= 0 Then
	        //     BJValue = Application.Max(K - S, 0)
	        // Else
	    		if (t <= 0l) {
	    			return Math.abs(k - s);
	    		}
	    		if (-b > q || q == 0) {
	    			return BlackScholes(type, s, k, t, v, r, q);
	    		}
	    		// BJValue = BjerkPrice(CallPutFlag, K, S, r - b, -b, Dividend, sigma, t)
	    		price = BjerkStensInternal(type, k, s, t, v, -b, r - b, q);
	    } else if ("C".equals(type)) {
	    		if (t <= 0l) {
				return Math.abs(s - k);
			}
	    		if (b > r || q == 0) {
	    			return BlackScholes(type, s, k, t, v, r, q);
	    		}
	        //  BJValue = BjerkPrice(CallPutFlag, S, K, r, b, Dividend, sigma, t)
	    		price = BjerkStensInternal(type, s, k, t, v, b, r, q);
	    }
		return price;
	}
	
	private static double BjerkStensInternal(String type, double s, double k, double t, double v, double b, double r, double q) {
		// NEED TO CHECK: in the money American options are priced at intrinsic value
		if (k < s) {
			// return Math.abs(s - k);
		}
		// sig2 = sigma * sigma
		double v2 = v * v;
	    // beta = (0.5 - b / sig2) + ((b / sig2 - 0.5) ^ 2 + 2 * r / sig2) ^ 0.5 'use eqn 6
		double beta = (0.5 - (b / v2)) + Math.pow( ( Math.pow(( (b / v2) - 0.5), 2) + (2 * r / v2) ), 0.5);

	    // 'evaluate trigger price X as in eq (10)
	    // B_inf = (beta / (beta - 1)) * K
	    	double betainfinity = ( beta / (beta -1) ) * k;	
	    // B0 = Application.Max(K, (r / Dividend) * K)
	    	double betazero = Math.max(k, (r / q) * k);
	    	
	    // h = -(b * t + 2 * sigma * t ^ 0.5) * (K * K / ((B_inf - B0) * B0))
	    	double h = -((b * t) + (2 * v * Math.pow(t, 0.5))) * ( (k * k) / ((betainfinity - betazero) * betazero)); 

	    // X = B0 + (B_inf - B0) * (1 - Exp(h))
	    	double X = betazero + ((betainfinity - betazero) * (1 - Math.exp(h)));
	    	
	    	/*
	    	 * 	    If X < S Then
	        BjerkPrice = Abs(S - K)
	    Else
	        alpha = (X - K) * X ^ (-beta)
	        'calculate call price as in eq (4)
	        tmp1 = alpha * S ^ beta
	        tmp2 = alpha * BjerkPhi(S, t, beta, X, X, sigma, r, b)
	        tmp3 = BjerkPhi(S, t, 1, X, X, sigma, r, b)
	        tmp4 = BjerkPhi(S, t, 1, K, X, sigma, r, b)
	        tmp5 = K * BjerkPhi(S, t, 0, X, X, sigma, r, b)
	        tmp6 = K * BjerkPhi(S, t, 0, K, X, sigma, r, b)
	        BjerkPrice = tmp1 - tmp2 + tmp3 - tmp4 - tmp5 + tmp6
	    End If
	    	 */
	    	if (X < s) {
	    		// price equals intrinsic value
	    		return Math.abs(s - k);
	    	} else {
	    		double alpha = (X -k) * Math.pow(X, -beta);
	    		double tmp1 = alpha * Math.pow( s, beta);
	    		double tmp2 = alpha * BjerkStensPhi(s, t, beta, X, X, v, r, b);
	    		double tmp3 = BjerkStensPhi(s, t, 1, X, X, v, r, b);
	    		double tmp4 = BjerkStensPhi(s, t, 1, k, X, v, r, b);
	    		double tmp5 = k * BjerkStensPhi(s, t, 0, X, X, v, r, b);
	    		double tmp6 = k * BjerkStensPhi(s, t, 0, k, X, v, r, b);
	    		double bjprice = tmp1 - tmp2 + tmp3 - tmp4 - tmp5 + tmp6;
	    		return bjprice;
	    	}

	}
	
	private static double BjerkStensPhi(double s, double t, double gamma, double h, double X, double v, double r, double b) {
	    // sig2 = sigma * sigma
		double v2 = v * v;
	    // K = 2 * b / sig2 + (2 * gamma - 1) 'use eqn (9)
	    	double K = ((2 * b) / v2) + (2 * gamma - 1);
	    // lambda = -r + gamma * b + 0.5 * gamma * (gamma - 1) * sig2 'use eqn (8)
	    	double lambda = -r + (gamma * b) + ((0.5 * gamma * (gamma - 1)) * v2);
	    // tmp1 = (Log(S / h) + (b + (gamma - 0.5) * sig2) * t) / (sigma * t ^ 0.5)
	    	double tmp1 = ((Math.log(s / h)) + (b + (gamma - 0.5) * v2) * t) / (v * Math.pow(t, 0.5));
	    // tmp2 = (Log(X ^ 2 / (S * h)) + (b + (gamma - 0.5) * sig2) * t) / (sigma * t ^ 0.5)
	    	double tmp2 = (Math.log((X * X) / (s * h)) + (b + (gamma - 0.5) * v2) * t) / (v * Math.pow(t, 0.5)); 
	    // tmp2 = tmp1 + 2 * Log(X / S) / (sigma * t ^ 0.5)
	    	tmp2 = tmp1 + 2 * Math.log(X / s) / (v * Math.pow(t, 0.5));
	    	
	    // BjerkPhi = Exp(lambda * t) * S ^ gamma * (Application.NormSDist(-tmp1) - ((X / S) ^ K) * Application.NormSDist(-tmp2))
	    	double psphi = Math.exp(lambda * t) * Math.pow(s, gamma) * (cdf(-tmp1) - (Math.pow(X / s, K)) * cdf(-tmp2));
		
		return psphi;

	}

}
