package lc1.stats;

import java.io.PrintWriter;

import lc1.util.Constants;
import pal.math.MultivariateFunction;
import pal.math.OrthogonalHints;
import cern.colt.matrix.DoubleMatrix2D;

public class OrthogonalProbabilityDistribution  implements lc1.stats.ProbabilityDistribution2, MultivariateFunction, Comparable{

	public lc1.stats.ProbabilityDistribution distx;
	public lc1.stats.ProbabilityDistribution disty;
	
	public OrthogonalProbabilityDistribution(lc1.stats.ProbabilityDistribution distx2,
			lc1.stats.ProbabilityDistribution disty2) {
		this.distx = distx2;
		this.disty = disty2;
		this.cuty = disty2 instanceof CutNormal;
	}
	
	
	public void setPriors(ProbabilityDistribution2 pr2, int type, boolean x
	){
	 if(x) this.distx.setPriors(((OrthogonalProbabilityDistribution)pr2).distx, type);	
	 else this.disty.setPriors(((OrthogonalProbabilityDistribution)pr2).disty, type);	
	}
final boolean cuty;

public void print(PrintWriter pw) {
	pw.print("x:\t");this.distx.print(pw);
	pw.println();
	pw.print("y:\t");this.disty.print(pw);
	
}

	
	public void addCount(double obj_index, double obj_indexy, double value) {
		if(value > Constants.countThresh2()){
			
		disty.addCount(obj_indexy, value);
		if(cuty ){
			double quantileL = value*((CutNormal)disty).quantileL;
			double quantileG = value*((CutNormal)disty).quantileG;
			double quantileRem = value*((CutNormal)disty).quantileRem;
			
			if(quantileL>Constants.countThresh2()){
				distx.addCount(obj_index, quantileL);
			}
			if(quantileG>Constants.countThresh2()){
				distx.addCount(obj_index, quantileG);
			}
			if(quantileRem>Constants.countThresh2()){
			distx.addCount(obj_index,quantileRem);
			}
		}
		else{
			distx.addCount(obj_index, value);
		}
		}
		
	}
	
	
	public int numObs() {
		return this.distx.numObs();
	}
	
	public int fill(DoubleMatrix2D x, DoubleMatrix2D y,DoubleMatrix2D yB, int numObs,double[] noCop,  double pseudo){
		
		
			int res1=   distx.fill(x, y, numObs, noCop, pseudo);
		///note values may not be in same order!
			if(yB!=null){
				int res2 =  disty.fill(x,yB,numObs, noCop,   pseudo);
				if(res2!=res1 && res2>=0){
					throw new RuntimeException(" !! "+res2+" "+res1+" ");//+((CutNormal)disty).maxToAdd.size()+" "+((CutNormal)disty).minToAdd.size());
				}
			}
		//	if(((TrainableNormal)distx).obsx.size()>0){
		//		System.err.println("h");
		//	}
			return res1;
	}
	
	 public int fillVariance(DoubleMatrix2D y, DoubleMatrix2D yb,
				DoubleMatrix2D covar, int numObs, double pseudo){
		 int res1 = distx.fillVariance(y,numObs, pseudo);
		 if(yb!=null){
			 int res2 = disty.fillVariance(yb,numObs, pseudo);
			 if(res2!=res1 && res2>=0) throw new RuntimeException("!!");
		 }
		return  res1;
	 }
	 
	 public void addCount(Double r2, Double b, double val,
				SimpleExtendedDistribution1 mixe1,
				ProbabilityDistribution2 probDistG) {
		 OrthogonalProbabilityDistribution dist1 = (OrthogonalProbabilityDistribution)probDistG;
		this.distx.addCount(r2,  val, mixe1, dist1.distx);
		this.disty.addCount( b, val, mixe1, dist1.disty);
	 }
	 
	 public void addCount(Double r2, Double b, double val,
				SimpleExtendedDistribution1 mixeR,SimpleExtendedDistribution1 mixeB,
				OrthogonalProbabilityDistribution dist1) {
	
		this.distx.addCount(r2,  val, mixeR, dist1.distx);
		this.disty.addCount( b, val, mixeB, dist1.disty);
	 }
	 public void variance( int type, double[] sum){
    	if(type==0) this.distx.variance(sum);
    	else disty.variance(sum);
    }
	
	
	public void setParam(int type,int i,double rho){
		if(type==0 || type==1) {
			if(i==0) {
				distx.setParam(type,rho);
			}
			else if(i==1){
			disty.setParam(type,rho);
			}
		}
	 }
	
	 
	 
	 
	 
	/*
	public void addCounts(ProbabilityDistribution2 probabilityDistribution) {
		distx.addCounts(((SkewNormal2)probabilityDistribution).distx);
		disty.addCounts(((SkewNormal2)probabilityDistribution).disty);
		
	}*/

	
	public ProbabilityDistribution2 clone(double u) {
		return new OrthogonalProbabilityDistribution(distx.clone(u), disty.clone(u));
	}

	
	
	public ProbabilityDistribution2 clone() {
		return new OrthogonalProbabilityDistribution(distx.clone(), disty.clone());
	}

	/*
	public double[] getCount(double[] angle) {
		// TODO Auto-generated method stub
		return null;
	}*/

	
	public void getIntervalX(double[] in, double[] res_inR) {
		distx.getInterval(in, res_inR);
		
	}

	
	public void getIntervalY(double[] in, double[] res_inR) {
		disty.getInterval(in, res_inR);
		
	}

	/*
	public double[] getMean() {
		return new double[] {
				distx.getMean(),disty.getMean()
				
		};
	}
	*/

	
	public int getParamIndex() {
		return distx.getParamIndex()+disty.getParamIndex();
	}

	/*
	public double getParamValue(int n1) {
	     if(n1<distx.getNumArguments()) return distx.getParamValue(n1);
	     else return disty.getParamValue(n1 - distx.getNumArguments());
	}*/

	
	public String id() {
		return distx.id()+" "+disty.id();
	}

	
	public void initialise() {
	distx.initialise();
	disty.initialise();	
		
	}

	
	public void maximise(double d1, double e1, double f1, double d2, double e2,
			double f2, double g) {
	//if(false){
		if(Constants.trainR())distx.maximise(d1, e1, f1);
	    if(Constants.trainB()) disty.maximise(d2, e2, f2);
	//}
		
	} 
	
	public String toString(){
		return distx.toString()+"x"+disty.toString();
	}

	
	public String name() {
		return distx.name()+disty.name();
	}

	/*
	public double prior() {
		throw new RuntimeException("!!");
	}*/

	
	public double probability(double x, double y) {
		double res;
		//if(Math.abs(y-0.7159317)<1e-4){
	//		System.err.println("res" );
	//	}
		if(Constants.suppressR()){
			int index =0;
			 res =  disty.probability(y) *  1.0/(Constants.maxR(index)- Constants.minR(index));//
			//	if(true) throw new RuntimeException("assuming index is 0");
		}
		else if(Constants.suppressB()){
			int index =0;
			//int index = ((PseudoDistribution)distx).data_index;
			res = distx.probability(x) * 1/(Constants.maxB(index)- Constants.minB(index));
			//if(true) throw new RuntimeException("assuming index is 0");
		}
		else if(Constants.rweight()<1.0){
			res =Math.pow( distx.probability(x),Constants.rweight()) * disty.probability(y);
		}
		else{
			res = distx.probability(x) * disty.probability(y);
				//Math.pow(distx.probability(x),0.02)*Math.pow(disty.probability(y),1.0);
		}
		
		return res;
	}

	
	public void recalcName() {
		distx.recalcName();
		disty.recalcName();
		
	}

/*	
	public double sample() {
	throw new RuntimeException("!!");
	}

	
	public void setParamValue(int n1, double val) {
		throw new RuntimeException("!!");
		
	}

	
	public void setParamsAsAverageOf(ProbabilityDistribution2[] tmp) {
		throw new RuntimeException("!!");
		
	}*/

	
	public void transfer(double pseudoC) {
		distx.transfer(pseudoC);
		disty.transfer(pseudoC);
		
	}

	/*
	public void transfercounts(EmissionState innerState, int phen_index, int i) {
		throw new RuntimeException("!!");
		
	}*/

	
	public void updateParamIndex() {
		distx.updateParamIndex();
		disty.updateParamIndex();
		
	}

	
	public int compareTo(Object o) {
		throw new RuntimeException("!!");
	}

	
	public double evaluate(double[] argument) {
		throw new RuntimeException("!!");
	}

	
	public double getLowerBound(int n) {
		throw new RuntimeException("!!");
		//return 0;
	}

	
	public int getNumArguments() {
		throw new RuntimeException("!!");
		//return 0;
	}

	
	public OrthogonalHints getOrthogonalHints() {
		throw new RuntimeException("!!");
	//	return null;
	}

	
	public double getUpperBound(int n) {
		throw new RuntimeException("!!");
		//return 0;
	}

	
	public void getInterval(double[] input, DoubleMatrix2D cov, double[] mean) {
	
		cov.setQuick(0, 0, Math.pow((this.distx).scale(), 2));
		cov.setQuick(1, 1, Math.pow((this.disty).scale(), 2));
		
		cov.setQuick(0, 1, 0);
		cov.setQuick(1, 0, 0);
		mean[0] = (this.distx).getMean();
		mean[1] =  (this.disty).getMean();
	}
    
	
	public ProbabilityDistribution2 clone(double u,
			SimpleExtendedDistribution1 dist1) {
		return new OrthogonalProbabilityDistribution(distx.clone(u, dist1),
				disty.clone(u, dist1));
		
	}
	
	public ProbabilityDistribution2 clone(double u, SimpleExtendedDistribution1 sx, SimpleExtendedDistribution1 sy){
		return new OrthogonalProbabilityDistribution(distx.clone(u, sx),
				disty.clone(u, sy));
	}
	
	public double probability(double x, double y, int mixComponent) {
		// TODO Auto-generated method stub
		//if(true) throw new RuntimeException("!!");
		double res;
		if(Constants.suppressR()){
			 res =  disty.probability(y, mixComponent) * 1;
			
		}
		else if(Constants.suppressB()){
			res = distx.probability(x,mixComponent);
//			int index = ((PseudoDistribution)distx).data_index;
//			res = distx.probability(x, mixComponent)* 1.0/(Constants.maxR(index)- Constants.minR(index));// * 1/(Constants.maxB()- Constants.minB());
		}
		else{
			res = distx.probability(x, mixComponent) * disty.probability(y, mixComponent);
				//Math.pow(distx.probability(x),0.02)*Math.pow(disty.probability(y),1.0);
		}
		
		return res;
	}
	
	public void setMinMax(double minR, double maxR, double min, double max) {
		distx.setMinMax(minR, maxR);
		disty.setMinMax(min,max);
		
	}
	
	public void setToExclude() {
		throw new RuntimeException("!!");
		
	}


	
}