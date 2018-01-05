package com.github.jasonmerecki.optpricing.black;

import static org.junit.Assert.*;

import org.junit.Test;

// http://www.fintools.com/resources/online-calculators/options-calcs/options-calculator/

public class BlackLikeTest {

	@Test
	public void testBlackCall1() {
		// Result       /  Online calculator
		// ---------------------------------------------
		// 20.037       / https://www.mystockoptions.com/black-scholes.cfm
		// 20.2961      / https://www.erieri.com/blackscholes
		// 20.2961667   / (excel spreadsheet)
		// 20.2961      / http://www.fintools.com/resources/online-calculators/options-calcs/options-calculator/
		
		double s = 1177.62d;
		double k = 1195.00d;
		double t = 0.084931506849315d; // date 12/19/2017, expiration 1/19/2018, 31 days
		double v = 0.20d;
		double r = 0.0135d;
		double q = 0.0d;
		double bsprice = BlackLike.priceBlackScholes("C", s, k, t, v, r, q);
		System.out.println("testBlackCall1 bsprice=" + bsprice);
		assertEquals(20.29616303951127d, bsprice, 0.00000000000d);
	}

	@Test
	public void testBlackPut1() {
		// Result       /  Online calculator
		// ---------------------------------------------
		// n/a          / https://www.mystockoptions.com/black-scholes.cfm
		// 0.2708       / https://www.erieri.com/blackscholes
		// ?????        / (excel spreadsheet)
		// 0,2708       / http://www.fintools.com/resources/online-calculators/options-calcs/options-calculator/
		
		double s = 214.76d;
		double k = 190.00d;
		double t = 0.084931506849315d; // date 12/19/2017, expiration 1/19/2018, 31 days
		double v = 0.25d;
		double r = 0.0135d;
		double q = 0.0d;
		double bsprice = BlackLike.priceBlackScholes("P", s, k, t, v, r, q);
		System.out.println("testBlackPut1 bsprice=" + bsprice);
		assertEquals(0.2707906395245452d, bsprice, 0.00000000000d);
	}
	
	@Test
	public void testBjerkStensCall1() {
		// Result       /  Online calculator
		// ---------------------------------------------
		// 19.082612   / (excel spreadsheet)

		double s = 1177.62d;
		double k = 1195.00d;
		double t = 0.084931506849315d; // date 12/19/2017, expiration 1/19/2018, 31 days
		double v = 0.20d;
		double r = 0.0135d;
		double q = 0.03d;
		double bsprice = BlackLike.priceBjerkStens("C", s, k, t, v, r, q);
		System.out.println("testBjerkStensCall1 bsprice=" + bsprice);
		assertEquals(19.082618995152643d, bsprice, 0.00000000000d);
	}
	
	@Test
	public void testBjerkStensPut1() {
		// Result       /  Online calculator
		// ---------------------------------------------
		// 22.0387792   / (excel spreadsheet)

		double s = 1177.62d;
		double k = 1165.00d;
		double t = 0.084931506849315d; // date 12/19/2017, expiration 1/19/2018, 31 days
		double v = 0.20d;
		double r = 0.0135d;
		double q = 0.03d;
		double bsprice = BlackLike.priceBjerkStens("P", s, k, t, v, r, q);
		System.out.println("testBjerkStensPut1 bsprice=" + bsprice);
		assertEquals(22.03875264497185d, bsprice, 0.00000000000d);
	}
	
}
