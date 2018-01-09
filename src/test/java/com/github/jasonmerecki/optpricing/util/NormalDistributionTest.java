package com.github.jasonmerecki.optpricing.util;

import static org.junit.Assert.*;

import org.apache.commons.math3.special.Erf;
import org.junit.Test;

public class NormalDistributionTest {
	
	@Test 
    public void testPdfStdMean() {
		NormalDistribution nd = new NormalDistribution();
        double x = 0.0d;
        double pdf = nd.pdf(x);
        double expectedPdf = 0.3989422804014327;
        double ulp = Math.ulp(expectedPdf); 
        System.out.println("for x=" + x + " the pdf=" + pdf + " expected=" + expectedPdf + " within " + ulp);
        assertEquals(expectedPdf, pdf, ulp);
	}
	
    @Test 
    public void testPdfStdNegativePositiveSame() {
        NormalDistribution nd = new NormalDistribution();
        double x = 0.56d;
        double pdf = nd.pdf(x);
        double expectedPdf = 0.34104578863035256d;
        double ulp = Math.ulp(expectedPdf); // see Math class javadoc
        System.out.println("for x=" + x + " the pdf=" + pdf + " expected=" + expectedPdf + " within " + ulp);
        assertEquals(expectedPdf, pdf, ulp);
        x = -0.56d;
        pdf = nd.pdf(x);
        System.out.println("for x=" + x + " the pdf=" + pdf + " expected=" + expectedPdf + " within " + ulp);
        assertEquals(expectedPdf, pdf, ulp);
    }
    
    @Test 
    public void testPdfStdTwoStdDev() {
    		NormalDistribution nd = new NormalDistribution();
        double x = 2.0d;
        double pdf = nd.pdf(x);
        double expectedPdf = 0.05399096651318806;
        double ulp = Math.ulp(expectedPdf); 
        System.out.println("for x=" + x + " the pdf=" + pdf + " expected=" + expectedPdf + " within " + ulp);
        assertEquals(expectedPdf, pdf, ulp);
	}
    
    
    @Test 
    public void testPdfMaxAtMean() {
    		NormalDistribution nd = new NormalDistribution();
    		double pdf = 0.0d;
    		double maxX = 0.0d;
    		for (double x = -8.0d; x <= 8.0d; x += 0.01d) {
    	        double testPdf = nd.pdf(x);
    	        // un-comment this line to observe the imprecision of the double type
    	        // System.out.println("testing x=" + x + " the testPdf=" + testPdf);
    	        if (testPdf > pdf) {
    	        		pdf = testPdf;
    	        		maxX = x;
    	        }
    		}
    		double expectedMaxX = 0.0d;
        double expectedPdf = 0.3989422804014327;
        double ulp = Math.ulp(expectedPdf); 
        System.out.println("for maxX=" + maxX + " the pdf=" + pdf + " expected=" + expectedPdf + " within " + ulp);
        assertEquals(expectedPdf, pdf, ulp);
        assertEquals(expectedMaxX, maxX, 0.000000000001d);
	}
    
    @Test
    public void testPdfVsApacheMath() {
    		NormalDistribution nd = new NormalDistribution();
        double x = 0.56d;
        double pdf = nd.pdf(x);
        org.apache.commons.math3.distribution.NormalDistribution and = new org.apache.commons.math3.distribution.NormalDistribution();
        double expectedPdf = and.density(x);
        double ulp = Math.ulp(expectedPdf); 
        System.out.println("for x=" + x + " the pdf=" + pdf + " expected(from Apache)=" + expectedPdf + " within " + ulp);
        assertEquals(expectedPdf, pdf, ulp);
    }
    
    
    @Test
    public void testPdfVsApacheMathNonstd() {
    		NormalDistribution nd = new NormalDistribution(34.3d, 3.22d);
        double x = 37.5d;
        double pdf = nd.pdf(x);
        org.apache.commons.math3.distribution.NormalDistribution and = new org.apache.commons.math3.distribution.NormalDistribution(34.3d, 3.22d);
        double expectedPdf = and.density(x);
        double ulp = Math.ulp(expectedPdf); 
        System.out.println("for mean/stdev 34.3d, 3.22d x=" + x + " the pdf=" + pdf + " expected(from Apache)=" + expectedPdf + " within " + ulp);
        assertEquals(expectedPdf, pdf, ulp);
    }
    
    @Test
    public void testCdf() {
		NormalDistribution nd = new NormalDistribution();
        double x = 0.56d;
        double cdf = nd.cdf(x);
        double expectedCdf = 0.7122603051006894d;
        double ulp = Math.ulp(expectedCdf); 
        System.out.println("for x=" + x + " the cdf=" + cdf + " expected=" + expectedCdf + " within " + ulp);
        assertEquals(expectedCdf, cdf, ulp);
    }
    
    @Test
    public void testCdf2() {
		NormalDistribution nd = new NormalDistribution(34.3d, 3.22d);
        double x = 37.8d;
        double cdf = nd.cdf(x);
        double expectedCdf = 0.8614719786451529d;
        double ulp = Math.ulp(expectedCdf); 
        System.out.println("for mean/stdev 34.3d, 3.22d, x=" + x + " the cdf=" + cdf + " expected=" + expectedCdf + " within " + ulp);
        assertEquals(expectedCdf, cdf, ulp);
    }
    
    @Test
    public void testCdfVsApache() {
		NormalDistribution nd = new NormalDistribution();
        double x = 0.56d;
        double cdf = nd.cdf(x);
        org.apache.commons.math3.distribution.NormalDistribution and = new org.apache.commons.math3.distribution.NormalDistribution();
        double expectedCdf = and.cumulativeProbability(x);
        double ulp = 1.0e-7; 
        System.out.println("for x=" + x + " the cdf=" + cdf + " expected(from Apache)=" + expectedCdf + " within " + ulp);
        assertEquals(expectedCdf, cdf, ulp);
    }
    
    @Test
    public void testCdfVsApache2() {
		NormalDistribution nd = new NormalDistribution(34.3d, 3.22d);
        double x = 37.8d;
        double cdf = nd.cdf(x);
        org.apache.commons.math3.distribution.NormalDistribution and = new org.apache.commons.math3.distribution.NormalDistribution(34.3d, 3.22d);
        double expectedCdf = and.cumulativeProbability(x);
        double ulp = 1.0e-7; 
        System.out.println("for mean/stdev 34.3d, 3.22d, x=" + x + " the cdf=" + cdf + " expected(from Apache)=" + expectedCdf + " within " + ulp);
        assertEquals(expectedCdf, cdf, ulp);
    }
    
    @Test
    public void testStdNormalCdf() {
        double x = 0.56d;
        double cdf = NormalDistribution.StandardNormal.cdf(x);
        double expectedCdf = 0.7122603051006894d;
        double ulp = Math.ulp(expectedCdf); 
        System.out.println("stdnormal, for x=" + x + " the cdf=" + cdf + " expected=" + expectedCdf + " within " + ulp);
        assertEquals(expectedCdf, cdf, ulp);
    }
    
    @Test
    public void testStdNormalCdfConvertShortcut() {
    		NormalDistribution nd = new NormalDistribution(34.3d, 3.22d);
    		double x = 37.5d;
    		// convert
    		double xconv = (x - 34.3d) / 3.22d;
        double cdfConvert = NormalDistribution.StandardNormal.cdf(xconv);
        double expectedCdf = nd.cdf(x);
        double ulp = Math.ulp(expectedCdf); 
        System.out.println("stdnormal, for x=" + x + " the cdfConvert=" + cdfConvert + " expected=" + expectedCdf + " within " + ulp);
        assertEquals(expectedCdf, cdfConvert, ulp);
    }
    
    @Test
    public void testStdNormalPdf() {
        double x = 0.56d;
        double pdf = NormalDistribution.StandardNormal.pdf(x);
        double expectedPdf = 0.34104578863035256d;
        double ulp = Math.ulp(expectedPdf); 
        System.out.println("stdnormal, for x=" + x + " the pdf=" + pdf + " expected=" + expectedPdf + " within " + ulp);
        assertEquals(expectedPdf, pdf, ulp);
    }
    
    @Test
    public void testErrFunction() {
    		double x = 0.56;
        double erf = NormalDistribution.errorFunction(x);
        double expectedErf = 0.5716157766617889d;
        double ulp = 1.0e-3; 
        System.out.println("erf, for x=" + 1 + " the erf=" + erf + " expectedErf=" + expectedErf + " within " + ulp);
        assertEquals(expectedErf, erf, ulp);
    }
    
    @Test
    public void testErrFunctionVsApache() {
    		double x = 0.56;
        double erf = NormalDistribution.errorFunction(x);
        double expectedErf = Erf.erf(0.56);
        double ulp = 1.0e-3; 
        System.out.println("erf, for x=" + 1 + " the erf=" + erf + " expectedErf=" + expectedErf + " within " + ulp);
        assertEquals(expectedErf, erf, ulp);
    }
    
}
