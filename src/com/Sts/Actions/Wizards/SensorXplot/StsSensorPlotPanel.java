package com.Sts.Actions.Wizards.SensorXplot;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import info.monitorenter.gui.chart.*;
import info.monitorenter.gui.chart.controls.*;
import info.monitorenter.gui.chart.pointpainters.*;
import info.monitorenter.gui.chart.traces.*;
import info.monitorenter.gui.chart.traces.painters.*;
import info.monitorenter.gui.chart.views.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.text.*;


/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSensorPlotPanel extends StsJPanel implements ChangeListener
{
    private StsSensorXplotWizard wizard;
    private StsSensorPlot wizardStep;
 
    private StsModel model = null;
    
    StsGroupBox axisBox = new StsGroupBox("");    
    
    StsGroupBox xyAxisBox = new StsGroupBox("Define X-Y Axis");
    StsGroupBox innerBox = new StsGroupBox();
    StsComboBoxFieldBean xAttributeBean = new StsComboBoxFieldBean(); 
    StsComboBoxFieldBean yAttributeBean = new StsComboBoxFieldBean();
    
    StsGroupBox colorBox = new StsGroupBox("Color");
    StsColorComboBoxFieldBean colorComboBox = null;

    StsButton clearBtn, clearCurrentBtn, newChartBtn, deleteChartBtn;
    StsBooleanFieldBean invertBean = new StsBooleanFieldBean();
    
	StsGroupBox messageBox = new StsGroupBox("Messages");
    StsGroupBox innerBox2 = new StsGroupBox();
	JTextPane msgBean = new JTextPane();
	StsButton exportBtn;
	StsButton clusterBtn;
	
    transient XplotChart[] charts = null;
    
    transient StsTimeCurve[] curves = null;
    transient StsColor stsColor = new StsColor(StsColor.RED);
    transient int colorIdx = 0;
    transient boolean inclusive = true;
    
	transient JTabbedPane tabbedPanels = new JTabbedPane();
	
    public StsSensorPlotPanel(StsSensorXplotWizard wizard, StsSensorPlot wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }
    
    public void buildPanel()
    {   
    	tabbedPanels.addChangeListener(this);
    	deleteChartBtn = new StsButton("Delete Level", "Delete current cross-plot level", this, "deleteChart");            	    	
    	newChartBtn = new StsButton("Add Level", "Add another cross-plot level", this, "addChart");            	
    	xyAxisBox.gbc.fill = gbc.HORIZONTAL;
    	innerBox.addToRow(newChartBtn);
        innerBox.addEndRow(deleteChartBtn);
        xyAxisBox.addEndRow(innerBox);

        xyAxisBox.gbc.fill = gbc.HORIZONTAL;
        xyAxisBox.gbc.anchor = gbc.WEST; 
    	xAttributeBean.initialize(this, "xAttribute", "X:", curves);
    	yAttributeBean.initialize(this, "yAttribute", "Y:", curves);
    	xyAxisBox.addEndRow(xAttributeBean);
    	xyAxisBox.addEndRow(yAttributeBean);
    	
        gbc.fill = gbc.HORIZONTAL;
        
        clearBtn = new StsButton("Clear All", "Clear all polygons and reset 3D view", this, "clearAll");        
        clearCurrentBtn = new StsButton("Clear Current", "Clear current color polygon and reset 3D view", this, "clearCurrent");
        invertBean.initialize(this, "inverted", "Include Events in Polygons", false);
        
        colorComboBox = new StsColorComboBoxFieldBean(this, "polygonColor", "Color:", StsColor.colors32);
        colorBox.gbc.anchor = gbc.WEST;
        colorBox.gbc.fill = gbc.HORIZONTAL;
        colorBox.gbc.gridwidth = 2;
        colorBox.addEndRow(colorComboBox);
        colorBox.addEndRow(invertBean);
        colorBox.gbc.fill = gbc.NONE;
        colorBox.gbc.gridwidth = 1;
        colorBox.addToRow(clearBtn);
        colorBox.addEndRow(clearCurrentBtn);
        
        exportBtn = new StsButton("Export View", "Export the sensor events currently in view.", wizard, "exportView");    	        
        clusterBtn = new StsButton("Save Clusters", "Save the current cluster values as a an attribute within the respective sensors", wizard, "saveClusters");        
    }
    
    public void stateChanged(ChangeEvent ev)
    {
    	xAttributeBean.doSetValueObject(charts[getChartIdx()].getXAttribute());
    	yAttributeBean.doSetValueObject(charts[getChartIdx()].getYAttribute());
    	revalidate();
    }
    
    public void addChart()
    {
    	buildChart(charts.length);
    	tabbedPanels.setSelectedIndex(charts.length-1);
    }
    public void deleteChart(int idx)
    {    
    	if(idx == 0)
    	{
    		new StsMessage(wizard.frame, StsMessage.WARNING, "Cannot delete the last level. A minimum of one level is required.");
    		return;
    	}
    	charts[idx] = null;
    	charts = (XplotChart[])StsMath.arrayDeleteElement(charts, idx);
    	tabbedPanels.remove(idx);
    	tabbedPanels.setSelectedIndex(charts.length-1);
    	updateView();   	
    }
    
    public void deleteChart()
    {
    	deleteChart(tabbedPanels.getSelectedIndex());
    } 
    
    public void revalidateNext()
    {
    	if(charts.length > getChartIdx()+1)
    		charts[getChartIdx()+1].rebuildTraces();
    }
    
    public boolean updateView()
    {
    	Main.logUsageTimer();
    	// Initialize sensor clusters
    	wizard.getPrimeSensor().setClustering(true);
    	int[] clusters = new int[wizard.getPrimeSensor().getNumValues()];
    	for(int i=0; i<clusters.length; i++)
    		clusters[i] = -1;
    		
    	// One the values in the last chart are valid.
    	int lastChartIdx = charts.length-1;
    	if(charts[lastChartIdx].valIndices == null)
        {
    	    if(lastChartIdx > 0)
    		    lastChartIdx--;
    	    else
    	    {
    		    wizard.getPrimeSensor().setClustering(false);
    		    wizard.getPrimeSensor().setClusters(null);
    		    return true;
    	    }
        }
    	
        int[] validIndices = charts[lastChartIdx].valIndices;
        int[] chartCluster = charts[lastChartIdx].clusters;
    	for(int i=0; i<validIndices.length; i++)
    	{
    		if(chartCluster[i] != -1)
    			clusters[validIndices[i]] = chartCluster[i];
    	}
    	
    	// Set sensor cluster array
    	wizard.getPrimeSensor().setClusters(clusters);
    	wizard.getModel().viewObjectRepaint(this, wizard.getPrimeSensor());
    	return true;
    }
    
    public boolean verifyEditOkay()
    {
    	int lastChartIdx = charts.length-1;
    	if(lastChartIdx != getChartIdx())
    	{
    		if(!StsYesNoDialog.questionValue(wizard.frame, "Changing polygons at higher level will result in lower levels being deleted.\nDo you want to continue?"))
    			return false;
    	}
    	return true;
    }
    
    public void deleteChartsAfterCurrent()
    {
    	int lastChartIdx = charts.length-1;
    	if(getChartIdx() >= lastChartIdx)
    		return;
    	for(int i=getChartIdx()+1; i<=lastChartIdx; i++)
    		deleteChart(i);
    }
    
    public void destroyChart()
    {
    	for(int i=0; i<charts.length; i++)
    	{
    		charts[i].destroyChart();
    		remove(charts[i]);
    	}
    }

    public void setInverted(boolean val)
    {
    	inclusive = val;
    	clearAll(false);
    	charts[getChartIdx()].reCompute();
    }
    
    public boolean getInverted() { return inclusive; }
    public void initialize()
    {
    	if(charts != null) return;
        buildChart(0);
        StsObject[] curves = wizard.getPrimeSensor().getPropertyCurves();
        curves = (StsObject[])StsMath.arrayAddElement(curves, wizard.getPrimeSensor().getRelativeTimeCurve());
    	
    	xAttributeBean.setListItems(curves);
    	yAttributeBean.setListItems(curves);
    	xAttributeBean.setSelectedItem(wizard.getPrimeSensor().getTimeCurve(StsLogVector.X));
    	yAttributeBean.setSelectedItem(wizard.getPrimeSensor().getTimeCurve(StsLogVector.Y));
    	
    	messageBox.gbc.fill = gbc.BOTH;
    	messageBox.add(msgBean);
    	messageBox.gbc.fill = gbc.HORIZONTAL;    	
    	innerBox2.addToRow(exportBtn);
    	innerBox2.addEndRow(clusterBtn);
        messageBox.addEndRow(innerBox2);
    	
        axisBox.addEndRow(xyAxisBox);
        axisBox.addEndRow(colorBox);
        axisBox.gbc.fill = gbc.BOTH;
        axisBox.addEndRow(messageBox);
        
        gbc.fill = gbc.VERTICAL;
        gbc.weightx = 0.0;
        gbc.anchor = gbc.WEST;
        addToRow(axisBox);
        
        gbc.weightx = 1.0;
        gbc.fill = gbc.BOTH;
        addEndRow(tabbedPanels);
        
        wizard.rebuild(); 
        
    	colorComboBox.setSelectedIndex(0);
    	stsColor = colorComboBox.getStsColor();
    	
        xAttributeBean.setBeanObject(this);
        yAttributeBean.setBeanObject(this);    	
    }
    
    public void setPolygonColor(StsColor color)
    {
    	stsColor = color;
    	colorIdx = colorComboBox.getSelectedIndex();
    }
    public int getChartIdx() 
    { 
    	if(tabbedPanels.getSelectedIndex() == -1)
    		return 0;
    	else
    		return tabbedPanels.getSelectedIndex(); 
    }
    
    public StsColor getPolygonColor()
    {
    	return stsColor;   	
    }
    
    public StsObject getXAttribute() { return charts[getChartIdx()].getXAttribute(); }
    public StsObject getYAttribute() { return charts[getChartIdx()].getYAttribute(); }
    public void setXAttribute(StsObject curve) { charts[getChartIdx()].setXAttribute(curve); }
    public void setYAttribute(StsObject curve) { charts[getChartIdx()].setYAttribute(curve); }    
    
    public void updateMessage(String msg)
    {
    	msgBean.setText(msg);
    }
    
    public void clearClusters()
    {
    	StsSensor sensor = wizard.getPrimeSensor();
    	sensor.setClustering(false);
    	sensor.setClusters(null);
    	wizard.getModel().viewObjectRepaint(this, sensor);
    }
    
    public void clearCurrent()
    {   	
    	charts[getChartIdx()].clearCurrent();
    	updateView();
    } 
    
    public void clearAll()
    {
    	clearAll(true);   	
    } 
    
    public void clearAll(boolean clearChart)
    {
    	if(clearChart) 
    		charts[getChartIdx()].clearAll();
    	wizard.getPrimeSensor().setClusters(null);
    	for(int i=charts.length-1; i>-1; i--)
    		charts[i].reCompute();
    }
    
    public void buildChart(int idx)
    {
        XplotChart newChart = new XplotChart(this, wizard.getPrimeSensor());
    	charts = (XplotChart[])StsMath.arrayAddElement(charts, newChart);
    	
        new LayoutFactory.BasicPropertyAdaptSupport(this.getRootPane(), charts[idx]);
        LayoutFactory lf = LayoutFactory.getInstance();
        lf.setShowTraceZindexMenu(false);
        lf.setShowTraceNameMenu(false);
        lf.setShowAxisXRangePolicyMenu(false);
        lf.setShowTraceNameMenu(true);

        charts[idx].setEnablePopup(false);
        charts[idx].setMinPaintLatency(100);
        
        charts[idx].setBorder(BorderFactory.createEtchedBorder());
    	tabbedPanels.add("Level #" + (idx+1), charts[idx]);        
    }
    
    public StsColor getClusterStsColor(int index)
    {
        return StsColor.colors32[index];
    }
}

class stsTracePainterDisc extends TracePainterDisc
{
	private stsPointPainterDisc m_pointPainter;	
	public stsTracePainterDisc(int size)
	{
		this.m_pointPainter = new stsPointPainterDisc(size);
    }	
}

class stsPointPainterDisc extends PointPainterDisc
{	
	public stsPointPainterDisc(int size)
	{
		super(size);
	}
	
	public void paintPoint(final int absoluteX, final int absoluteY, final int nextX,
		      final int nextY, final Graphics2D g, final TracePoint2D original) 
	{
        g.drawRect(absoluteX - m_halfDiscSize, (absoluteY - m_halfDiscSize), m_discSize, m_discSize);
    }
}

class XplotChart extends Chart2D implements MouseListener, MouseMotionListener 
{
    transient DecimalFormat labelFormat = new DecimalFormat("###0.0");
    transient StsSensorPlotPanel panel = null;
    transient Polygon[] m_area = new Polygon[StsColor.colors32.length];
    transient Point2D startPt;
    
    int[] valIndices = null;
    public int[] clusters = null;
    ChartPanel chartPanel = null;
    StsSensor sensor = null;
    transient ITrace2D trace = null;
    transient ITrace2D htrace = null;
    transient int previousIdx = 0;
    StsTimeCurve xAttribute = null;
    StsTimeCurve yAttribute = null;
    
    boolean debug = false;
    
    public XplotChart(StsSensorPlotPanel panel, StsSensor sensor) 
    {
		super();
		xAttribute = (StsTimeCurve)sensor.getTimeCurve(StsLogVector.X);
		yAttribute = (StsTimeCurve)sensor.getTimeCurve(StsLogVector.Y);
		this.sensor = sensor;
		this.panel = panel;		
		buildTraces();		
		chartPanel = new ChartPanel(this);	
		addMouseListener(this);
		addMouseMotionListener(this);
	}
    
    public void rebuildTraces()
    {
    	/* Convert current polys to real coordinates */
    	double[][][] polys = new double[m_area.length][][];
    	for(int i=0; i<m_area.length; i++)
    	{
    		if(m_area[i]!= null)
    		{
        		polys[i] = new double[m_area[i].npoints][2];    			
    			int[] xVals = m_area[i].xpoints;
    			int[] yVals = m_area[i].ypoints;
    			for(int j=0; j<m_area[i].npoints; j++)
    			{
    				polys[i][j][0] = getAxisX().translatePxToValue(xVals[j]);
    				polys[i][j][1] = getAxisY().translatePxToValue(yVals[j]);
    			}
    		}
    	}     	
    	removeTrace(trace);
    	buildTraces();
    	/* Convert polys to new screen coordinates */
    	for(int i=0; i<m_area.length; i++)
    	{
    		if(m_area[i]!= null)
    		{
    			int nPts = m_area[i].npoints;
    			for(int j=0; j<nPts; j++)
    			{
    				m_area[i] = new Polygon();
    				m_area[i].addPoint(getAxisX().translateValueToPx(polys[i][j][0])
    						, getAxisY().translateValueToPx(polys[i][j][1]));
    			}
        		updateSensor(m_area[i],i);    			
    		}
    	}
    	setRequestedRepaint(true);
    }
    
    private void buildTraces()
    {
        Color color = Color.BLACK;
        trace = createDataset();
        trace.setColor(color);
        addTrace(trace);
        trace.setName("Level #" + panel.getChartIdx());
    } 
    
    private Trace2DLtd createDataset()
    {
        Trace2DLtd trace = new Trace2DLtd(sensor.getNumValues());
        
        ITracePainter painter = new TracePainterDisc(4);
        trace.setTracePainter(painter);
        int[] sClusters = sensor.getClusters();
        
        float[]	xVals = sensor.getPropertyAsFloats(xAttribute);
        float[]	yVals = sensor.getPropertyAsFloats(yAttribute);
        valIndices = new int[xVals.length];
        clusters = new int[xVals.length];
        int cnt = 0;
        for (int i = 0; i < xVals.length; i++)
        {
            if((xVals[i] == StsParameters.nullValue) || (yVals[i] == StsParameters.nullValue))
                continue;
        	if(sClusters != null)
        	{
            	if(sClusters[i] != -1)
            	{
            		trace.addPoint(xVals[i], yVals[i]);
            		valIndices[cnt] = i;
            		clusters[cnt++] = -1;
            	}
        	}
        	else
        	{
        		trace.addPoint(xVals[i], yVals[i]);
        		valIndices[cnt] = i;
        		clusters[cnt++] = -1;
        	}
        }
        valIndices = (int[])StsMath.trimArray(valIndices, cnt);
        clusters = (int[])StsMath.trimArray(clusters, cnt);
        return trace;
    }
    
    public StsObject getXAttribute() { return xAttribute; }
    public StsObject getYAttribute() { return yAttribute; }
    public void setXAttribute(StsObject curve)
    {
    	if(xAttribute == curve) 
    		return;
    	
    	if(!panel.verifyEditOkay())
    		return;
    	else
    		panel.deleteChartsAfterCurrent();
    	
    	xAttribute = (StsTimeCurve)curve;
        
    	clearAll();
    	valIndices = null;
    	clusters = null;    	
    	panel.updateView();    	
    	removeTrace(trace);
    	buildTraces();
        setRequestedRepaint(true);
    }
    public void setYAttribute(StsObject curve)
    {
    	if(yAttribute == curve) 
    		return;
    	
    	if(!panel.verifyEditOkay())
    		return;
    	else
    		panel.deleteChartsAfterCurrent();
    	
    	yAttribute = (StsTimeCurve)curve;
    	
    	clearAll();
    	valIndices = null;
    	clusters = null;
    	panel.updateView();    	
    	removeTrace(trace);
    	buildTraces();    	
        setRequestedRepaint(true);    	
    }
    
    public void highlightPoint()
    {
    	if(sensor.getPicked() == -1) return;
    	
    	if(htrace != null)
    	{
    		removeTrace(htrace);
    		htrace = null;
    	}
    	
    	htrace = new Trace2DLtd(3);
    	ITracePainter painter = new stsTracePainterDisc(8);
    	htrace.setTracePainter(painter);
    	
        float[]	xVals = sensor.getPropertyAsFloats(xAttribute);
        float[]	yVals = sensor.getPropertyAsFloats(yAttribute);
        
        htrace.setColor(Color.RED);
        htrace.addPoint(trace.getMaxX(), trace.getMaxY());
        htrace.addPoint(trace.getMinX(), trace.getMinY());
        htrace.addPoint(xVals[sensor.getPicked()], yVals[sensor.getPicked()]);
        htrace.setName("Selected");
        
        addTrace(htrace);
        setCurrentTrace(0);
        return;
    }
    
    public void destroyChart()
    {
        trace = null;
    }
    public void updateSensor(int cIdx)
    {
    	if(m_area[cIdx]!= null)
    		updateSensor(m_area[cIdx], cIdx);
    }   
    
    public void updateSensor(Polygon poly, int cIdx)
    {
    	if(poly.npoints < 3) 
    		return;
    	
    	if(!panel.verifyEditOkay())
    		return;
    	else
    		panel.deleteChartsAfterCurrent();
    	
        float[]	xVals = sensor.getPropertyAsFloats(xAttribute);
        float[]	yVals = sensor.getPropertyAsFloats(yAttribute);
    	for(int i=0; i<valIndices.length; i++)
    	{	
    		int xpx = getAxisX().translateValueToPx(xVals[valIndices[i]]);
    		int ypx = getAxisY().translateValueToPx(yVals[valIndices[i]]);
    		if(poly.contains(xpx, ypx))
    		{
    			if(panel.inclusive)
    				clusters[i] = cIdx;
    			else
    				clusters[i] = -1;
    		}
    	}
    	if(debug)
    	{
    		for(int i=0; i<clusters.length; i++)
    			System.out.println("Cluster[" + i + "]=" + clusters[i]);
    	}    	
    	panel.updateView();
    }    
    	// Mouse Listeners
     public void mouseExited(MouseEvent e) { }
     public void mouseEntered(MouseEvent e) 
     { 
    	 highlightPoint();
     }
     public void mouseReleased(MouseEvent e)
     {
         double mouseY = (double) getAxisY().translatePxToValue(e.getY());
         double mouseX = (double) getAxisX().translatePxToValue(e.getX());
         String msg = xAttribute.getName() + "=" + labelFormat.format(mouseX) + "\n" +
            	yAttribute.getName() + "=" + labelFormat.format(mouseY);
         panel.updateMessage(msg);
         
         return;         
     }
     public void mouseClicked(MouseEvent e) {  } 
     public void clearAll()
     {
    	 for(int i=0; i<m_area.length; i++)
    		 m_area[i] = null;
    	 for(int i=0; i<clusters.length; i++)
    			 clusters[i] = -1;
         if(panel.getChartIdx() == 0)
         	sensor.setClusters(null);
         setRequestedRepaint(true);
         panel.updateView();
     }
     
     public void clearCurrent()
     {
    	 m_area[panel.colorIdx] = null;
    	 for(int i=0; i<clusters.length; i++)
    	 {
    		 if(clusters[i] == panel.colorIdx)
    			 clusters[i] = -1;
    	 }
         setRequestedRepaint(true);   		 
     }
     
     public void reCompute()
     {
         for(int i=0; i<m_area.length; i++)
         {
         	if(m_area[i] != null)
         		updateSensor(m_area[i], i);
         }  	 
     }
     public void mousePressed(MouseEvent e)
     {  
         if(m_area[panel.colorIdx] == null)
        	 m_area[panel.colorIdx] = new Polygon();
         m_area[panel.colorIdx].addPoint(e.getX(), e.getY());
         updateSensor(m_area[panel.colorIdx], panel.colorIdx);
         setRequestedRepaint(true);
     } 

     // Motion Listeners
     public void mouseDragged(MouseEvent e)
     {
        double mouseY = (double) getAxisY().translatePxToValue(e.getY());
        double mouseX = (double) getAxisX().translatePxToValue(e.getX());
        String msg = xAttribute.getName() + "=" + labelFormat.format(mouseX) + "\n" +
        	yAttribute.getName() + "=" + labelFormat.format(mouseY);
        panel.updateMessage(msg);       
     }
     public void mouseMoved(MouseEvent e) { } 
     public void paintComponent(final Graphics g) 
     {
        super.paintComponent(g);

        Graphics2D g2 = (Graphics2D) g;

        for(int i=0; i<m_area.length; i++)
        {
        	if(m_area[i] != null)
        	{
        		g2.setPaint(Color.BLACK);        		
        		g2.draw(m_area[i]);
        		g2.setPaint(StsColor.colors32[i].getColor());
        		g2.fill(m_area[i]);
        	}
        }
     }    
}

