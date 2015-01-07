package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Base class used to present two-dimensional data.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;
import info.monitorenter.gui.chart.*;
import info.monitorenter.gui.chart.axis.*;
import info.monitorenter.gui.chart.controls.*;
import info.monitorenter.gui.chart.labelformatters.*;
import info.monitorenter.gui.chart.traces.*;
import info.monitorenter.gui.chart.traces.painters.*;
import info.monitorenter.gui.chart.views.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import java.text.*;
import java.util.*;
import java.sql.Array;

public class StsViewTimeSeries extends JFrame implements MouseListener, MouseMotionListener, PropertyChangeListener, ComponentListener, StsSerializable, WindowListener
{
    static public final String viewClassname = "TimeSeries_View";
    static public byte SENSORS = 0;
    static public byte LIVEWELLS = 1;
    public byte plotType = SENSORS;
    public Object[] selectedObjects = null;
	public StsSensorClass sensorClass = null;
    public int currentSensor = 0;
    public Point windowLocation = null;
    public Dimension windowSize = null;
    public boolean[] visible = null;

    transient Chart2D chart = null;
    transient ChartPanel chartPanel = null;
    //transient JScrollPane chartScrollPane = new JScrollPane();
    transient ITrace2D[] traces = null;
    transient long maxTime = 0l;
    transient long minTime = StsParameters.largeLong;
    transient boolean axisChanged = true;

    transient float[][] attributes = null;
    transient int nAttributes = 0;
    transient String[] names = null;
    transient long[] time = null;
    transient boolean debug = false;
    transient long startMouseX = -1;
    transient long endMouseX = -1;
    transient float startMouseY = -1;
    transient float endMouseY = -1;
    transient boolean dragging = false;
    transient boolean axisSet = false;
	transient byte setting = BOTH;
    transient double maxY = -StsParameters.largeDouble;

    transient JPanel jPanel1 = new JPanel();
    transient JLabel rtLabel = null;

    transient GridBagLayout gridBagLayout2 = new GridBagLayout();
    transient GridBagLayout gridBagLayout1 = new GridBagLayout();

    transient StsJPanel panel = new StsJPanel();
    transient StsWin3dBase window = null;
    transient StsModel model = null;
    transient StsProject project = null;
    transient long processedTime = 0l;
    transient int axisType = -1;

    transient boolean setByMe = false;

    static final byte X = 0;
    static final byte Y = 1;
    static final byte BOTH = 2;
    static final byte INVALID = 3;
    static final byte NONE = -1;

    static final String[] invisibleAttributes = new String[] {"X", "Y", "Z", "DEPTH" };

    static public final String viewName = "TimeSeries_View";

    /**
     * Default constructor
     */
    public StsViewTimeSeries()
    {
    //    System.out.println("Default Time Series View Constructor.");
    }

    public StsViewTimeSeries(StsWin3dBase window)
    {
        super("Time Series Plot");
        // super(window, "Time Series Plot", false);
        try
        {
            this.window = window;
            this.model = window.getModel();
            this.project = model.getProject();
            currentSensor = 0;
            setAlwaysOnTop(true);
            StsSensorClass sensorClass = (StsSensorClass)model.getStsClass("com.Sts.DBTypes.StsSensor");
            StsLiveWellClass liveWellClass = (StsLiveWellClass)model.getStsClass("com.Sts.DBTypes.StsLiveWell");

            StsObject[] sensorObjects =  new StsObject[0];
            if(sensorClass != null)
                sensorObjects = sensorClass.getSensors();
            StsObject[] liveWellObjects = new StsObject[0];
            if(liveWellClass != null)
                liveWellObjects = liveWellClass.getObjectList();

            StsObject[] activeObjects = null;

            // Determine if plot is for wells or sensors if both are in project
            if((liveWellObjects.length > 0) && (sensorObjects.length > 0))
            {
                StsRadioButtonDialog dialog = new StsRadioButtonDialog(window,"Select Type", true, new String[] {"Sensors", "Live Wells"});
                dialog.setVisible(true);
                String selected = dialog.getSelectedButton();
                if(selected.equalsIgnoreCase("Sensors"))
                {
                    plotType = SENSORS;
                    activeObjects = getActiveSensorObjects(sensorObjects);
                }
                else
                {
                    plotType = LIVEWELLS;
                    activeObjects = getActiveWellObjects(liveWellObjects);
                }
            }
            else if((liveWellObjects.length > 0) && (sensorObjects.length == 0))
            {
                plotType = LIVEWELLS;
                activeObjects = getActiveWellObjects(liveWellObjects);
            }
            else
            {
                plotType = SENSORS;
                activeObjects = getActiveSensorObjects(sensorObjects);
            }

            if(this instanceof StsViewRealtimeSeries)
            {
                if(plotType == LIVEWELLS)
                {
                    new StsMessage(window, StsMessage.WARNING, "Real-time time series are not available for wells.");
                    return;
                }
                selectedObjects = new Object[1];
                selectedObjects[0] = StsListSelectionDialog.getSingleSelectFromListDialog(window, "Selection", null, sensorObjects);
            }
            else
            {
                while(!verifyAttributes(selectedObjects))
                {
                    selectedObjects = StsListSelectionDialog.getMultiSelectFromListDialog(window, "Selection", null, activeObjects);
                    if(selectedObjects == null)
                        break;
                }
            }
            if(selectedObjects == null)
            	return;
            if(selectedObjects[0] == null)
            	return;
            if(!buildDialog())
            {
                selectedObjects = null;
                return;
            }
            setSize(1000, 250);
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private boolean verifyAttributes(Object[] sensors)
    {
        if(sensors == null) return false;
        if(sensors.length == 1) return true;
        StsTimeCurve[] curves = ((StsSensor)sensors[0]).getPropertyCurvesExcludingXYZ();
        boolean isDynamic = ((StsSensor)sensors[0]).canBeDynamic();
        for(int i=1; i<sensors.length; i++)
        {
            if(((StsSensor)sensors[i]).canBeDynamic() != isDynamic)
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "Cannot mix static and dynamic sensors on the same plot....");
                return false;
            }
            if(!((StsSensor)sensors[i]).hasCurves(curves))
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "The same curves must exist in all selected sensors....");
                return false;
            }
        }
        return true;
    }

    private StsObject[] getActiveSensorObjects(StsObject[] objects)
    {
        StsObject[] newObjects = new StsObject[objects.length];
        int cnt = 0;
        for(int i=0; i<objects.length; i++)
        {
            if(((StsSensor)objects[i]).getNumValues() > 0)
                newObjects[cnt++] = objects[i];
        }
        return (StsObject[])StsMath.trimArray(newObjects, cnt);
    }

    private StsObject[] getActiveWellObjects(StsObject[] objects)
    {
        StsObject[] newObjects = new StsObject[objects.length];
        int cnt = 0;
        for(int i=0; i<objects.length; i++)
        {
            if(((StsLiveWell)objects[i]).getNPoints() > 0)
                newObjects[cnt++] = objects[i];
        }
        return (StsObject[])StsMath.trimArray(newObjects, cnt);
    }
    public boolean buildDialog()
    {
        try
        {
        	if(selectedObjects[0] == null) return false;
            this.setTitle("Time Series Plot - "  + selectedObjects[0].toString());

            if(!getAttributeData())
            {
                // Open empty time series view for realtime.
                new StsMessage(window, StsMessage.ERROR,  "No data found in selected object.\n\n\nIf you are starting a realtime series, \n" +
                        " start realtime using the toolbar controls\n and then launch plot again.");
                return false;
            }
            verifyRange();
            jbInit();
            initialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
            return false;
        }
        return true;
    }

    public void initializeTransients(StsWin3dBase window, StsModel model)
    {
        this.window = window;
        this.model = model;
        this.project = model.getProject();
        currentSensor = 0;
        //setAlwaysOnTop(true);
        start();
    }

    public boolean start()
    {
        try
        {
            buildDialog();
            restoreWindowState();
            setVisible(true);
            return true;
        }
        catch(Exception e)
        {
            e.printStackTrace();
            return false;
        }
    }

    public void initialize()
    {
        traces = new Trace2DLtd[nAttributes];

        chart = createChart();
        reconfigureView();

        chart.setTimeDateFormat(project.getTimeDateFormatString());
        chart.setMinPaintLatency(100);

        new LayoutFactory.BasicPropertyAdaptSupport(this.getRootPane(), chart);
        LayoutFactory lf = LayoutFactory.getInstance();
        lf.setShowTraceZindexMenu(false);
        lf.setShowTraceNameMenu(false);
        lf.setShowAxisXRangePolicyMenu(true);
        chart.setEnablePopup(false);
        chartPanel = new ChartPanel(chart);
        jPanel1.add(chartPanel,   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));

        chart.addMouseListener(this);
        chart.addMouseMotionListener(this);
        chart.addPropertyChangeListener(Chart2D.PROPERTY_PROPERTIES_PANEL, this);
        chart.setBackground(project.getTimeSeriesBkgdColor().getColor());
    	chart.setForeground(project.getTimeSeriesFrgdColor().getColor());
        chart.setGridColor(project.getTimeSeriesGridColor().getColor());
        this.addComponentListener(this);
        this.addWindowListener(this);
        axisChanged = true;
        reconfigureView();

        // Set current trac to first visible
        for(int i=0; i<visible.length; i++)
        {
            if(visible[i])
            {
                setCurrentTrace(i);
                break;
            }
        }
    }

    public void changePolicy(Chart2D chart1)
    {
        if(chart1 == null) return;
		Iterator traces = chart1.getTraces().iterator();
		while(traces.hasNext())
            ((Trace2DLtd) traces.next()).setMaxSize(getNumValues());

    	chart1.setRequestedRepaint(true);
        chart1.updateUI();
    }

    public int getNumValues()
    {
        int num = 0;
        for(int i=0; i<selectedObjects.length; i++)
        {
            if(plotType == SENSORS)
                num = num + ((StsSensor)selectedObjects[i]).getNumValues();
            else
                num = num + ((StsLiveWell)selectedObjects[i]).getNPoints();
        }
        return num;
    }
    public void setCurrentTraceToAmplitudeOrPressure()
    {
        for(int i=0; i<names.length; i++)
        {
            if((names[i].toLowerCase().contains("amp") || names[i].toLowerCase().contains("mag"))
                    || names[i].toLowerCase().contains("pres"))
            {
                setCurrentTrace(i);
                return;
            }
        }
    }

    public void destroyChart()
    {
        chartPanel = null;
        chart.destroy();
        chart = null;
        time = null;
        for(int i=0; i<traces.length; i++)
            traces[i] = null;
        traces = null;
        attributes = null;
        visible = null;
        names = null;
    }

    void jbInit() throws Exception
    {
        addWindowListener(this);
        panel.setLayout(gridBagLayout1);

        jPanel1.setLayout(gridBagLayout2);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());

        panel.add(jPanel1,   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));

        //chartScrollPane.getViewport().add(panel, null);
        //this.getContentPane().add(chartScrollPane);
        this.getContentPane().add(panel);
        this.addPropertyChangeListener(this);
    }

    private Chart2D createChart()
    {
        Chart2D chart1 = new Chart2D();
        buildTraces(chart1);
        chart1.setPaintLabels(true);
        chart1.setPaintCrosshair(true);
        return chart1;
    }

    public void buildTraces(Chart2D chart1)
    {
        boolean ampSet = false;
        ITracePainter barPainter = new TracePainterVerticalBar(2, chart1);
        Stroke thickLine = new BasicStroke(2);
        for(int i=0; i<names.length; i++)
        {
            Color color = StsColor.colors32[i%32].getColor();
            traces[i] = createDataset(names[i], attributes[i]);
            traces[i].addPropertyChangeListener(ITrace2D.PROPERTY_VISIBLE, this);
            traces[i].setColor(color);
            traces[i].setVisible(visible[i]);
            if((names[i].toLowerCase().contains("amp") || names[i].toLowerCase().contains("mag")) && !ampSet)
            {
                traces[i].setTracePainter(barPainter);
                ampSet = true;
            }
            else
                traces[i].setStroke(thickLine);
            
            //System.out.println("Setting number of samples to: " + getNumValues());
            ((Trace2DLtd)traces[i]).setMaxSize(getNumValues());
            chart1.addTrace(traces[i]);
            if(chart1.getCurrentTrace().getMaxY() != chart1.getCurrentTrace().getMinY())
                configureAxis(chart1);
        }
        //setCurrentTraceToAmplitude();
    	chart1.setRequestedRepaint(true);
        chart1.updateUI();
        changePolicy(chart1);
    }

    public void configureAxis(Chart2D chart1)
    {
        IAxis axis = chart1.getAxisX();
        axis.setPaintGrid(true);
        axis.setFormatter(new LabelFormatterDate(new SimpleDateFormat("kk:mm:ss")));

        axis = chart1.getAxisY();
        NumberFormat numFormat = NumberFormat.getNumberInstance();
        numFormat.setMaximumFractionDigits(2);
        numFormat.setMinimumFractionDigits(2);
        axis.setFormatter(new LabelFormatterNumber(numFormat));
        axis.setPaintGrid(true);
    }

    public void setVisibility(ITrace2D trace)
    {
    	for(int i=0; i<invisibleAttributes.length; i++)
    	{
    		if(trace.getName().equalsIgnoreCase(invisibleAttributes[i]))
    			trace.setVisible(false);
    	}
    }

    public void setCurrentTrace()
    {
    	if(chart == null) return;
        if(plotType == SENSORS)
        {
    	    for(int i=0; i<chart.getTraces().size(); i++)
    	    {
    		    if(((Trace2DLtd)chart.getTraces().toArray()[i]).isVisible())
                {
    			    chart.setCurrentTrace(i);
                }
    	    }
        }
    }

    public void setCurrentTrace(int idx)
    {
    	if(chart == null) return;
        visible[idx] = true;
        traces[idx].setVisible(true);
        chart.setCurrentTrace(idx);
        chart.updateUI();
    }

    public void propertyChange(final PropertyChangeEvent evt)
    {
        String property = evt.getPropertyName();
        if(visible == null) return;
        if (property.equals(ITrace2D.PROPERTY_VISIBLE))
        {
            int idx = getAttributeIndex(((ITrace2D)evt.getSource()).getLabel());
            visible[idx] = ((Boolean)evt.getNewValue()).booleanValue();
            if(visible[idx])
                setCurrentTrace(idx);      // Last activated trace
            else
                setCurrentTrace();         // Last live trace
            if(chart != null) configureAxis(chart);
        }
        if(chart != null) chart.propertyChange(evt);
    }

    public boolean reconfigureView()
    {
        byte newAxisType = StsSensor.LINEAR;
        if(plotType == SENSORS)
        {
            StsSensor sensor = (StsSensor)selectedObjects[0];
            newAxisType = sensor.getAxisType();
        }

        if(!axisChanged)
            return true;
        if(axisType != newAxisType)
        {
    	    switch(newAxisType)
    	    {
    	    case StsSensor.LOG10:
    		    chart.setAxisY(new AxisLog10());
    		    break;
    	    case StsSensor.LINEAR:
    		    chart.setAxisY(new AxisLinear());
    		    break;
    	    case StsSensor.LOGE:
    		    chart.setAxisY(new AxisLogE());
    		    break;
    	    }
            axisType = newAxisType;
            configureAxis(chart);            
        }
    	chart.setBackground(project.getTimeSeriesBkgdColor().getColor());
    	chart.setForeground(project.getTimeSeriesFrgdColor().getColor());
    	chart.setGridColor(project.getTimeSeriesGridColor().getColor());

        chart.getAxisX().setPaintGrid(project.getEnableTsGrid());
    	chart.getAxisY().setPaintGrid(project.getEnableTsGrid());

    	chart.setAreaCalculationIndex(0);
    	return true;
    }

    public int getAttributeIndex(String name)
    {
        for(int i=0; i<names.length; i++)
        {
            if(name.equalsIgnoreCase(names[i]))
                return i;
        }
        return -1;
    }
    
    public Trace2DLtd createDataset(String name, float[] attribute)
    {
        Trace2DLtd trace = new Trace2DLtd(500000);
        trace.setName(name);
        float min = StsMath.minExcludeNull(attribute);
        for (int i = 0; i < time.length; i++)
        {
        	if(maxTime < time[i]) maxTime = time[i];
        	if(minTime > time[i]) minTime = time[i];
            if((attribute[i] == model.getProject().getLogNull()) || (attribute[i] == StsParameters.nullValue))
                trace.addPoint(time[i], min);
            else
                trace.addPoint(time[i], attribute[i]);
        }
        if(trace.getMaxY() > maxY)
            maxY = trace.getMaxY();
        return trace;
    }

    public boolean initializeMovie()
    {
        return true;
    }

    public boolean getAttributeData()
    {
        if(plotType == SENSORS)
            return getSensorAttributeData();
        else
            return getLiveWellAttributeData();
    }

    public boolean getSensorAttributeData()
    {
        names = ((StsSensor)selectedObjects[0]).getTimeCurveList();
        nAttributes = names.length;
        if(nAttributes == 0) return false;
        attributes = new float[nAttributes][];
        int cnt = 0;
        if(visible == null)
        {
            visible = new boolean[nAttributes];
            for(int i=0; i<nAttributes; i++)
            {
                visible[i] = false;
                // Turn on amplitude
                if(names[i].toLowerCase().contains("amp") || names[i].toLowerCase().contains("mag"))
                {
                	visible[i] = true;
                    cnt++;
                }
                // Turn on obvious pump curve atts
                if(names[i].toLowerCase().contains("rate") || names[i].toLowerCase().contains("con") || names[i].toLowerCase().contains("prop") || names[i].toLowerCase().contains("pres"))
                {
                    visible[i] = true;
                    cnt++;
                }
            }
            if((nAttributes > 0) && (cnt == 0))
                visible[0] = true;  // Want something turned on
        }
        time = null;

        for(int j=0; j<selectedObjects.length; j++)
        {
            StsObjectRefList timeCurves = ((StsSensor)selectedObjects[j]).getTimeCurves();
            for(int i=0; i<timeCurves.getSize(); i++)
            {
                StsTimeCurve curve = (StsTimeCurve)timeCurves.getElement(i);
                if(i==0)
                {
                    //System.out.println("Adding " + curve.getNumValues() + " values from " + ((StsSensor)sensors[j]).getName());
                    time = StsMath.addPrimativeArrays(time, curve.getTimeVectorLongs());
                }
                attributes[i] = StsMath.addPrimativeArrays(attributes[i],curve.getValuesVectorFloats());
            }
        }
        return true;
    }

    public boolean getLiveWellAttributeData()
    {
        nAttributes = 2;   // Depth and Mdepth
        if(nAttributes == 0) return false;
        attributes = new float[nAttributes][];
        names = new String[] {"MDepth", "Depth"};        

        int cnt = 0;
        if(visible == null)
        {
            visible = new boolean[nAttributes];
            for(int i=0; i<nAttributes; i++)
                visible[i] = true;
        }
        time = null;

        for(int j=0; j<selectedObjects.length; j++)
        {
            StsLiveWell well = ((StsLiveWell)selectedObjects[j]);
            time = well.getTimeVector().getLongs();
            attributes[0] = StsMath.addPrimativeArrays(attributes[0],well.getMDepthFloats());
            attributes[1] = StsMath.addPrimativeArrays(attributes[1],well.getDepthFloats());
        }
        return true;
    }

    public void verifyRange()
    {
        if(plotType != SENSORS)
            visible[0] = true;
        else
        {
            // If max and min are the same set visible to false
            for(int k=0; k<attributes.length; k++)
            {
                float max = StsMath.max(attributes[k]);
                float min = StsMath.min(attributes[k]);
                if(min == max)
                    visible[k] = false;
            }

            // Make sure at least one curve is visible
            int cnt = 0;
            for(int i=0; i<visible.length; i++)
            {
                if(visible[i])
                {
                    cnt = 1;
                    break;
                }
            }
            if(cnt == 0) visible[0] = true;
        }
    }

    public void mouseExited(MouseEvent e) { }
    public void mouseEntered(MouseEvent e) { }
    public void mouseReleased(MouseEvent e)
    {
        if(e.getButton() == e.BUTTON3)
        {
        	// determine if a range is defined.
            if(endMouseX > startMouseX)
            	chart.zoomXRange(startMouseX, endMouseX);
            else
            	chart.unzoom();
            chart.setRequestedRepaint(true);
        	return;
        }   	
    	axisSet = false;
    	dragging = false;
        startMouseX = -1;
        if (chart != null)
        {
            double mouseY = (double) chart.getAxisY().translatePxToValue(e.getY());
            long mouseX = (long) chart.getAxisX().translatePxToValue(e.getX());
            byte setting = getActiveAxis(mouseX, mouseY);
            switch(setting)
            {
                case X:
                    endMouseX = mouseX;
                    endMouseX = checkXBounds(endMouseX);
                    setTime(endMouseX);
                    project.setProjectTime(endMouseX, true);
                    break;
                case Y:
                    break;
                case BOTH:
                    endMouseX = mouseX;
                    endMouseX = checkXBounds(endMouseX);
                    setTime(endMouseX);
                    project.setProjectTime(endMouseX, true);
                    break;
                case NONE:
                    break;
            }

        }
    }
    public void mouseClicked(MouseEvent e)
    {
    }
    public void mousePressed(MouseEvent e)
    {
    	if(chart != null)
        {
            double mouseY = (double) chart.getAxisY().translatePxToValue(e.getY());
            long mouseX = (long) chart.getAxisX().translatePxToValue(e.getX());
            byte setting = getActiveAxis(mouseX, mouseY);
            switch(setting)
            {
                case X:
                	if(endMouseX > startMouseX)
                	{
                		endMouseX = mouseX;
                	}
                	else
                	{
                		startMouseX = endMouseX = mouseX;
                		startMouseX = checkXBounds(startMouseX);
                		project.clearProjectTimeDuration();
                	}
                    setByMe = false;
                    setTime(startMouseX);
                    project.setProjectTime(endMouseX, true);
                    break;
                case Y:
                    break;
                case BOTH:
                    chart.clearXColorRange();
                    startMouseX = endMouseX = mouseX;
                    startMouseX = checkXBounds(startMouseX);
                    project.clearProjectTimeDuration();
                    setByMe = false;
                    setTime(startMouseX);
                    project.setProjectTime(endMouseX, true);
                    break;
                case INVALID:
                    StsMessageFiles.infoMessage("The min and max of the last selected attribute are equal, select a different attribute or turn off the last selected.");
                    break;
                case NONE:
                    StsMessageFiles.infoMessage("Set ranges by clicking on X or Y axis.");
                    break;
            }
        }
    }

    private byte getActiveAxis(long x, double y)
    {
        long xMin = (long)chart.getMinX();
        long xMax = (long)chart.getMaxX();
        float yMin = (float)chart.getMinY();
        float yMax = (float)chart.getMaxY();

        if(yMin == yMax)
            return INVALID;
        else if(y <= yMin && x >= xMin)
            return X;
        else if(x <= xMin && y >= yMin)
            return Y;
        else if(x < xMax && x > xMin &&  y < yMax && y > yMin)
            return BOTH;
        else
            return NONE;
    }

    public void mouseDragged(MouseEvent e)
    {
        long duration = 0;
        if(chart != null)
        {
            double mouseY = (double) chart.getAxisY().translatePxToValue(e.getY());
            long mouseX = (long) chart.getAxisX().translatePxToValue(e.getX());
            if(!axisSet)
            	setting = getActiveAxis(mouseX, mouseY);
            switch(setting)
            {
                case X:
                    //chart.clearXColorRange();
                    if(endMouseX > startMouseX)
                    {
                    	duration = endMouseX - startMouseX;
                        endMouseX = mouseX;
                        startMouseX = endMouseX - duration;
                        chart.setXColorRange(startMouseX, endMouseX);
                    }
                    else
                    {
                        if(startMouseX == -1) return;
                    	startMouseX = endMouseX = mouseX;
                    }
                    
                    setTime(endMouseX);
                    project.setProjectTime(endMouseX, true);
                    setByMe = true;
                    axisSet = true;
                    break;
                case Y:
//                    System.out.println("Y");
//                    if(startMouseY == -1) return;
//                    dragging = true;
//                    endMouseY = mouseY;
//                    endMouseY = checkYBounds(endMouseY);
//                    if(endMouseY < startMouseY)
//                    {
//                        model.clearTimeViews();
//                        chart.setYColorRange(startMouseY, endMouseY);
//                        // Set Sensor Range
//                        int idx = getAttributeIndex(chart.getCurrentTrace().getLabel());
//                        sensor.getSensorClass().setValueRangeAndSensor(sensor, idx, endMouseY, startMouseY);
//                    }
//                    setYValue(endMouseY);
//                    setByMe = true;
                    break;
                case BOTH:
                    if(startMouseX != -1)
                    {
                        dragging = true;
                        endMouseX = mouseX;
                        endMouseX = checkXBounds(endMouseX);
                        if(endMouseX > startMouseX)
                        {
                            chart.setXColorRange(startMouseX, endMouseX);
                            project.setProjectTimeDuration(startMouseX, endMouseX);
                        }
                        project.setProjectTime(endMouseX, true);
                        setTime(endMouseX);
                        setByMe = true;
                    }
//                    if(startMouseY != -1)
//                    {
//                        endMouseY = mouseY;
//                        endMouseY = checkYBounds(endMouseY);
//                        if(endMouseY < startMouseY)
//                        {
//                            model.clearTimeViews();
//                            chart.setYColorRange(startMouseY, endMouseY);
//                            // Set Sensor Range
//                            int idx = getAttributeIndex(chart.getCurrentTrace().getLabel());
//                            sensor.getSensorClass().setValueRangeAndSensor(sensor, idx, endMouseY, startMouseY);
//                        }
//                        setByMe = true;
//                        setYValue(endMouseY);
//                    }
                    break;
                case NONE:
                    break;
            }
            chart.setRequestedRepaint(true);
        }
    }

    public void clearView()
    {
        chart.clearYColorRange();
    }

    public void resetView(StsProject project)
    {
    	chart.setTimeDateFormat(project.getTimeDateFormatString());
    }

    public long checkXBounds(long value)
    {
//        System.out.println("Check X=" + value + " Max=" + chart.getMaxX() + " Min=" + chart.getMinX());
        if(value > chart.getMaxX())
            value = (long) chart.getMaxX();
        if(value < chart.getMinX())
            value = (long)chart.getMinX();
        return value;
    }

    public float checkYBounds(float value)
    {
//        System.out.println("Check Y=" + value + " Max=" + chart.getMaxY() + " Min=" + chart.getMinY());
        if(value >= chart.getMaxY())
            value = (float) chart.getMaxY();
        if(value <= chart.getMinY())
            value = (float)chart.getMinY();
        return value;
    }

    public void mouseMoved(MouseEvent e)
    {

    }

    public boolean viewObjectChanged(Object source, Object object)
    {
    	if(object instanceof StsTimeActionToolbar)
    	{
            //changePolicy(chart);
    		return true;
    	}
    	else if(object instanceof StsSensor)
        {
        	if(isPlottedSensor((StsSensor)object))
        	{
                reconfigureView();
        	}
        }
    	else if(object instanceof StsSensorClass)
        {
            reconfigureView();
        }
    	return false;
    }

    public boolean isPlottedSensor(StsSensor object)
    {
        for(int i=0; i<selectedObjects.length; i++)
        {
            if(selectedObjects[i] == object)
                return true;
        }
        return false;
    }

    public boolean viewObjectRepaint(Object source, Object object)
    {
        if(object instanceof StsProject)
        {
            StsProject project = (StsProject)object;
            setTime(project.getProjectTime());
        }
        return true;
    }

    public void setTime(long v_)
    {
        final long v = v_;
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                long value = v;
                if(debug) System.out.println("Set time to - " + value);
                if((value < 0) || (chart == null))
                    return;
                long duration = value - project.getProjectTimeDuration();
                if((duration > chart.getMaxX()) && (setByMe))
                {
                    setByMe = false;
                    project.clearProjectTimeDuration();
                    chart.clearXColorRange();
                }
                value = checkXBounds(value);
                chart.setCrosshairXValue(value);
                if(project.getProjectTimeDuration() != 0)
                    chart.setXColorRange(duration, value);
                else
                    chart.clearXColorRange();
            }
        };
        StsToolkit.runWaitOnEventThread(runnable);
    }

    public void setYValue(float value)
    {
        if((value < 0) || (chart == null)) return;
        value = checkYBounds(value);
        chart.setCrosshairYValue(value);
    }

    public void windowOpened(WindowEvent e) { }
    public void windowClosing(WindowEvent e)
    {
        project.clearProjectTimeDuration();
        if(plotType == SENSORS)
            ((StsSensor)selectedObjects[0]).getSensorClass().setValueRangeAndSensor(null, 0, 0, 0);
        model.deleteAuxTimeWindow(this);
    }
    public void windowClosed(WindowEvent e) {}
    public void windowIconified(WindowEvent e) {}
    public void windowDeiconified(WindowEvent e) {}
    public void windowActivated(WindowEvent e) {}
    public void windowDeactivated(WindowEvent e) {}
    public void componentShown(ComponentEvent e) {}
    public void componentResized(ComponentEvent e)
    {
        windowSize = new Dimension(getWidth(), getHeight());
    }
    public void componentMoved(ComponentEvent e)
    {
        windowLocation = new Point(getX(), getY());
    }
    public void componentHidden(ComponentEvent e) {}
    public void restoreWindowState()
    {
    	if(chart == null) // Chart has not been restored yet.
    		return;
    	if(windowSize == null)
        {
            if(!chart.getTraces().isEmpty())
                setSize(1000,200);
            else
    		    return;   // Rogue window is getting persisted.
        }
        else
        {
            setSize(windowSize);
            setLocation(windowLocation);
        }
        if(visible != null)
        {
            int cnt = 0;
            Iterator it = chart.getTraces().iterator();
            ITrace2D trace;
            while (it.hasNext())
            {
                trace = (ITrace2D) it.next();
                trace.setVisible(visible[cnt]);
                cnt++;
            }
        }
        reconfigureView();
        setEnabled(true);
    }

    public Frame getFrame()
    {
        return (Frame)this.getOwner();
    }

    public Class[] getViewClasses() { return new Class[] { StsViewTimeSeries.class }; }
}
