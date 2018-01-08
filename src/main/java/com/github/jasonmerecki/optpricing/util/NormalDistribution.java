package com.github.jasonmerecki.optpricing.util;

// references: 
// http://mathworld.wolfram.com/NormalDistribution.html
// http://mathworld.wolfram.com/RegularizedGammaFunction.html
// http://mathworld.wolfram.com/Erf.html
// https://introcs.cs.princeton.edu/java/21function/ErrorFunction.java.html

public class NormalDistribution {
	private static final double TOO_MANY_DEVIATIONS = 8d; // Apache commons Math uses 40.
	private final double mean;
	private final double standardDeviation;
	private final double stdDevSqrtTwoPi;
	private final double twoSigmaSquared;
	private final double sqrt2 = Math.sqrt(2.0d);
	
	public NormalDistribution(double mean, double standardDeviation) {
		this.mean = mean;
		this.standardDeviation = standardDeviation;
		this.stdDevSqrtTwoPi = standardDeviation * (Math.sqrt(2d * Math.PI));
		this.twoSigmaSquared = standardDeviation == 1d ? 2d : 2d * Math.pow(standardDeviation, 2d);
	}
	
	public NormalDistribution() {
		this(0, 1);
	}
	
	public double pdf(double x) {
		double expon;
		if (mean == 0d) {
			expon = -(x*x) / twoSigmaSquared; 
		} else {
			expon = -(Math.pow((x-mean),2d)) / twoSigmaSquared;
		}
		double probDist = Math.exp(expon) / stdDevSqrtTwoPi;
		return probDist;
	}
	
	public double cdf(double x) {
		final double distFromMean = x - mean;
        if (Math.abs(distFromMean) > TOO_MANY_DEVIATIONS * standardDeviation) {
        		return x < mean ? 0.0d : 1.0d;
        }
        	
        double errf = errorFunction(distFromMean / (standardDeviation * sqrt2));
		double cdf = 0.5d * (1.0 + errf);
		return cdf;
	}
	
	public static double errorFunction(double z) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(z));

        // use Horner's method
        double ans = 1 - t * Math.exp( -z*z   -   1.26551223 +
                                            t * ( 1.00002368 +
                                            t * ( 0.37409196 + 
                                            t * ( 0.09678418 + 
                                            t * (-0.18628806 + 
                                            t * ( 0.27886807 + 
                                            t * (-1.13520398 + 
                                            t * ( 1.48851587 + 
                                            t * (-0.82215223 + 
                                            t * ( 0.17087277))))))))));
        if (z >= 0) return  ans;
        else        return -ans;
    }

	public double getMean() {
		return mean;
	}

	public double getStandardDeviation() {
		return standardDeviation;
	}
	
	// static calculations for the standard normal - helps reduce object creation
	public static class StandardNormal {
		private static final NormalDistribution STD_NORM = new NormalDistribution();
		// cannot instantiate
		private StandardNormal() {}
		public static double getMean() { return STD_NORM.getMean(); }
		public static double getStandardDeviation() { return STD_NORM.getStandardDeviation(); }
		public static double pdf(double x) {
			return STD_NORM.pdf(x);
		}
		public static double cdf(double x) {
			return STD_NORM.cdf(x);
		}
		public static double cdf(double x, double mean, double standardDeviation) {
			return STD_NORM.cdf( ((x - mean) / standardDeviation) );
		}
	}
	
}
