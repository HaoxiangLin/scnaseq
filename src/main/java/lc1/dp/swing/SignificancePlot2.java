package lc1.dp.swing;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.logging.Logger;

import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;

import lc1.dp.core.BaumWelchTrainer;
import lc1.dp.core.ProbMultivariate;
import lc1.dp.core.Sampler;
import lc1.dp.data.collection.ArmitageCalculator;
import lc1.dp.data.collection.AssociationCalculator;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.LinearRegressionCalculator;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.model.CompoundMarkovModel;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.stats.Mixture;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.SkewNormal;
import lc1.stats.StateDistribution;
import lc1.util.Constants;

import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.AbstractVectorGraphicsIO;
import org.freehep.graphicsio.ImageConstants;
import org.freehep.graphicsio.ImageGraphics2D;
import org.freehep.graphicsio.emf.EMFGraphics2D;
import org.freehep.graphicsio.pdf.PDFGraphics2D;
import org.freehep.graphicsio.svg.SVGGraphics2D;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.event.AxisChangeEvent;
import org.jfree.chart.event.AxisChangeListener;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.CategoryItemRenderer;
import org.jfree.chart.renderer.category.DefaultCategoryItemRenderer;
import org.jfree.chart.renderer.xy.AbstractXYItemRenderer;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYBubbleRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYZDataset;


public class SignificancePlot2 {
    
  
   public static File makeOrDelete(File outdir, String st){
       File fi = new File(outdir, st);
       if(!fi.exists()) fi.mkdir();
       else{
       File[] f = fi.listFiles();
       for(int i=0; i<f.length; i++){
           f[i].delete();
       }
       }
       return fi;
   }
 
 static int divsize = 1;
private String[] names2;
//private ColorAdapter[] ca_rstates;

public static void init(List<String> snpid, int len){
	if(snp_alias!=null) return;
	
		String[] snpsToPlot = new String[snpid.size()];
     	probesToPlot = snpsToPlot;
     	snp_alias = new int[probesToPlot.length];
     	for(int i=0; i<snp_alias.length; i++){
     		snp_alias[i]=  i;
     		probesToPlot[i] = snpid.get(i);
     	}
    // 	pvalueToP = new double[len][snpsToPlot.length];
     //	minId = new int[len];
    start_ =0;
    end_ = snpsToPlot.length-1;
	 
}

  
   
    
    
    
    public static void reset(){
    	snp_alias = null;
    	probesToPlot = null;
    	//pvalueToP = null;
    }
    public  static int[] snp_alias;
    public static String[] probesToPlot;
   // public static double[][] pvalueToP;
    public static int start_, end_;
   
  


    public static final List<ProbabilityDistribution> extract(Map<String, List<ProbabilityDistribution>> m){
        List<ProbabilityDistribution> res = new ArrayList<ProbabilityDistribution>();
        for(Iterator<List<ProbabilityDistribution>> it = m.values().iterator(); it.hasNext();){
            res.add(it.next().get(0));
        }
        return res;
    }
    public static final Map<String, List<ProbabilityDistribution>> transform(Collection<ProbabilityDistribution> s2){
        Map<String, List<ProbabilityDistribution>> m = new HashMap<String, List<ProbabilityDistribution>>();
        
        for(Iterator<ProbabilityDistribution> it = s2.iterator(); it.hasNext();){
            ProbabilityDistribution pdist = it.next();
            String st = pdist.toString();
            List<ProbabilityDistribution > l = m.get(st);
            if(l==null){
                m.put(pdist.toString(), l = new ArrayList<ProbabilityDistribution>());
            }
            l.add(pdist);
        }
        return m;
    }
    
   static boolean rescale = false;
public static double plotThresh = 1e-3;

//publi

	// AbstractVectorGraphicsIO     gB,gR;
    final static int numPerPage = 10;
	
	public static int snp_alias(int i) {
		if(snp_alias==null) return i;
		int res = snp_alias[i];
		return res;
	}

    
    
}
