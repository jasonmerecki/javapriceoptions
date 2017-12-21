package com.github.jasonmerecki.optpricing.black;

import static org.junit.Assert.*;

import org.junit.Test;

// http://www.fintools.com/resources/online-calculators/options-calcs/options-calculator/

public class BlackLikeTest {

	@Test
	public void testBlackCall1() {
		// 26.65	 option price actual
		// 19.48% implied volatility
		
		// spot = 1,177.62
		// strike = 1195
		// January 19, 2018 expiry
		// dec 20, 2017 today
		// 12 + 19 = 31 days, = 0.084931506849315
		// volatility 0.20
		// rate = 0.0135
		double s = 1177.62d;
		double k = 1195.00d;
		double t = 0.084931506849315d;
		double v = 0.20d;
		double r = 0.0135d;
		double q = 0.0d;
		double bsprice = BlackLike.BlackScholes("C", s, k, t, v, r, q);
		System.out.println("dldldl bsprice=" + bsprice);
	}

}
