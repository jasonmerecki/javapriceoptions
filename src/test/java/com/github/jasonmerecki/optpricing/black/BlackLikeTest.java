package com.github.jasonmerecki.optpricing.black;

import static org.junit.Assert.*;

import org.junit.Test;

// http://www.fintools.com/resources/online-calculators/options-calcs/options-calculator/

public class BlackLikeTest {

	@Test
	public void testBlackCall1() {
		// 26.65	 option price actual
		// 19.48% implied volatility according to Yahoo! Finance
		
		// results of various calculators
		// http://www.fintools.com/resources/online-calculators/options-calcs/options-calculator/
		// fintools 19.8426 call
		// https://www.mystockoptions.com/black-scholes.cfm
		// mystockoptions 20.296
		// https://www.erieri.com/blackscholes
		// erieri.com 20.2961
		// Excel spreadsheet  20.2961667
		
		double s = 1177.62d;
		double k = 1195.00d;
		double t = 0.084931506849315d; // date 12/19/2017, expiration 1/19/2018, 31 days
		double v = 0.20d;
		double r = 0.0135d;
		double q = 0.0d;
		double bsprice = BlackLike.BlackScholes("C", s, k, t, v, r, q);
		System.out.println("dldldl bsprice=" + bsprice);
		assertEquals(19.75236950057615d, bsprice, 0.00000000000d);
	}

}
