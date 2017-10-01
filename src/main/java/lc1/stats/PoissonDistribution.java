package lc1.stats;

import java.io.PrintWriter;
import java.io.Serializable;

import lc1.dp.states.EmissionState;

import org.apache.commons.math.distribution.PoissonDistributionImpl;
import org.jfree.data.xy.XYSeries;

import pal.math.OrthogonalHints;
import cern.colt.matrix.DoubleMatrix2D;

public class PoissonDistribution implements ProbabilityDistribution, Serializable {

	private double rate;
	static PoissonDistributionImpl poisson;
	
	
	public PoissonDistribution(double rate)
	{
		this.rate=rate;
		poisson = new PoissonDistributionImpl(rate);
	}

	public PoissonDistribution(PoissonDistribution poissonDistribution, double rate) {
		this(rate);
	}

	
	public double probability(double x) {
		return poisson.probability(x);
	}

	
	public double getMean() {
		return poisson.getMean();
	}
	
	public void setCoverage(double depth) {
		poisson.setMean(depth+0.0001);
	}
	
	
	public void addCount(double objIndex, double value) {
		// TODO Auto-generated method stub
		
	}

	
	public void addCount(Double b, double val,
			SimpleExtendedDistribution1 mixe1, ProbabilityDistribution disty) {
		// TODO Auto-generated method stub
		
	}

	
	public void addCounts(ProbabilityDistribution probabilityDistribution) {
		// TODO Auto-generated method stub
		
	}

	
	public PoissonDistribution clone(){
		return new PoissonDistribution(rate);
	}
	
	
	public PoissonDistribution clone(double rate){
		return null;
	}

	
	public ProbabilityDistribution clone(double u,
			SimpleExtendedDistribution1 dist1) {
		// TODO Auto-generated method stub
		return null;
	}

	
	public int fill(DoubleMatrix2D x, DoubleMatrix2D y, int numObs,
			double[] noCop,  double pseudo) {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public int fillVariance(DoubleMatrix2D y, int numObs, double pseudo) {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public double[] getCount(double[] angle) {
		// TODO Auto-generated method stub
		return null;
	}

	
	public void getInterval(double[] in, double[] resInR) {
		// TODO Auto-generated method stub
		
	}



	
	public int getParamIndex() {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public double getParamValue(int n1) {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public String id() {
		// TODO Auto-generated method stub
		return null;
	}

	
	public void initialise() {
		// TODO Auto-generated method stub
		
	}

	
	public void maximise(double d, double e, double f) {
		// TODO Auto-generated method stub
		
	}

	
	public String name() {
		// TODO Auto-generated method stub
		return null;
	}

	
	public int numObs() {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public double plotObservations(String string, boolean b, XYSeries obs,
			boolean swtch) {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public void plotTheoretical(String string, boolean b, XYSeries theor) {
		// TODO Auto-generated method stub
		
	}

	
	public void print(PrintWriter pw) {
		// TODO Auto-generated method stub
		
	}

	
	public double prior() {
		// TODO Auto-generated method stub
		return 0;
	}


	
	public double probability(double x, int mixComponent) {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public void recalcName() {
		// TODO Auto-generated method stub
		
	}

	
	public double sample() {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public double scale() {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public void setMinMax(double min, double max) {
		// TODO Auto-generated method stub
		
	}

	
	public void setParam(int type, double rho) {
		// TODO Auto-generated method stub
		
	}

	
	public void setParamValue(int n1, double val) {
		// TODO Auto-generated method stub
		
	}

	
	public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
		// TODO Auto-generated method stub
		
	}

	
	public void setPriors(ProbabilityDistribution distx, int type) {
		// TODO Auto-generated method stub
		
	}

	
	public double sum() {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public void transfer(double pseudoC) {
		// TODO Auto-generated method stub
		
	}

	
	public void transfercounts(EmissionState innerState, int phenIndex, int i) {
		// TODO Auto-generated method stub
		
	}

	
	public void updateParamIndex() {
		// TODO Auto-generated method stub
		
	}

	
	public void variance(double[] sum) {
		// TODO Auto-generated method stub
		
	}

	
	public int compareTo(Object o) {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public double evaluate(double[] arg0) {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public double getLowerBound(int arg0) {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public int getNumArguments() {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public OrthogonalHints getOrthogonalHints() {
		// TODO Auto-generated method stub
		return null;
	}

	
	public double getUpperBound(int arg0) {
		// TODO Auto-generated method stub
		return 0;
	}


	

}
