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
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;
import info.monitorenter.gui.chart.*;
import info.monitorenter.gui.chart.controls.*;
import info.monitorenter.gui.chart.errorbars.ErrorBarPainter;
import info.monitorenter.gui.chart.errorbars.ErrorBarPolicyRelative;
import info.monitorenter.gui.chart.labelformatters.*;
import info.monitorenter.gui.chart.pointpainters.PointPainterDisc;
import info.monitorenter.gui.chart.pointpainters.PointPainterLine;
import info.monitorenter.gui.chart.traces.*;
import info.monitorenter.gui.chart.traces.painters.TracePainterDisc;
import info.monitorenter.gui.chart.traces.painters.TracePainterFill;
import info.monitorenter.gui.chart.traces.painters.TracePainterLine;
import info.monitorenter.gui.chart.traces.painters.TracePainterVerticalBar;
import info.monitorenter.gui.chart.views.*;
import info.monitorenter.util.Range;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import java.text.*;
import java.util.*;

public class StsViewTimeLine extends JDialog implements MouseListener, MouseMotionListener, PropertyChangeListener, ComponentListener, StsSerializable, WindowListener
{
    static public final String viewClassname = "TimeLine_View";
    public Point windowLocation = null;
    public Dimension windowSize = null;
    public boolean[] visible = null;

    transient Chart2D chart = null;
    transient ChartPanel chartPanel = null;
    transient ITrace2D[] traces = null;
    transient long maxTime = 0l;
    transient long minTime = StsParameters.largeLong;

    static byte STATICSENSOR = 0;
    static byte DYNAMICSENSOR = 1;
    static byte LIVEWELL = 2;
    static byte[] TIME_TYPES = new byte[] {STATICSENSOR, DYNAMICSENSOR, LIVEWELL};
    static String[] TIME_TYPE_STRINGS = new String[] {"StaticSensor", "DynamicSensor", "LiveWell"};

    transient float[][] attributes = null;
    transient long[][] times = null;
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

    transient JPanel jPanel1 = new JPanel();

    transient GridBagLayout gridBagLayout2 = new GridBagLayout();
    transient GridBagLayout gridBagLayout1 = new GridBagLayout();

    transient StsJPanel panel = new StsJPanel();
    transient JLabel timeLabel = new JLabel("Time:");
    transient StsWin3dBase window = null;
    transient StsModel model = null;
    transient StsProject project = null;
    transient long processedTime = 0l;

    transient boolean setByMe = false;

    static final byte X = 0;
    static final byte Y = 1;
    static final byte BOTH = 2;
    static final byte NONE = -1;

    static final String[] invisibleAttributes = new String[] {"X", "Y", "Z", "DEPTH" };

    static public final String viewName = "TimeLine_View";

    /**
     * Default constructor
     */
    public StsViewTimeLine()
    {
        System.out.println("Default Time Line View Constructor.");
    }

    public StsViewTimeLine(StsWin3dBase window)
    {
        super(window, "Time Line Plot", false);
        try
        {
            this.window = window;
            this.model = window.getModel();
            this.project = model.getProject();
            buildDialog();
            setSize(800, 250);
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }
    public StsViewTimeLine(StsWin3dBase window, boolean child)
    {
        super(window, "Time Line Plot", false);
        try
        {
            this.window = window;
            this.model = window.model;
            this.project = model.getProject();
            setSize(800, 250);
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void buildDialog()
    {
        try
        {
            if(!getAttributeData())
            {
                new StsMessage(window, StsMessage.WARNING, "Failed to load all time line data from database");
                return;
            }
            jbInit();
            initialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initializeTransients(StsWin3dBase window, StsModel model)
    {
        this.window = window;
        this.model = model;
        this.project = model.getProject();
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
        chart.setTimeDateFormat(project.getTimeDateFormatString());
        chart.setMinPaintLatency(100);

        new LayoutFactory.BasicPropertyAdaptSupport(this.getRootPane(), chart);
        LayoutFactory lf = LayoutFactory.getInstance();
        lf.setShowTraceZindexMenu(false);
        lf.setShowTraceNameMenu(false);
        lf.setShowAxisXRangePolicyMenu(false);
        chart.setEnablePopup(false);
        chartPanel = new ChartPanel(chart);
        jPanel1.add(chartPanel,   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));

        chart.addMouseListener(this);
        chart.addMouseMotionListener(this);
        chart.addPropertyChangeListener(Chart2D.PROPERTY_PROPERTIES_PANEL, this);
        this.addComponentListener(this);
        this.addWindowListener(this);
        reconfigureView();

        processedTime = (long)traces[0].getMaxX();        
    }

    public void changePolicy(Chart2D chart1)
    {
    	if(chart1 == null) return;
    	chart1.setRequestedRepaint(true);
        chart1.updateUI();
    }

    public void destroyChart()
    {
        time = null;
        chart = null;
        traces = null;
        attributes = null;
        repaint();
    }

    void jbInit() throws Exception
    {
        addWindowListener(this);
        panel.setLayout(gridBagLayout1);

        jPanel1.setLayout(gridBagLayout2);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());

        panel.add(jPanel1,   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
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

    private void buildTraces(Chart2D chart1)
    {
        Color color = null;
        ITracePainter linePainter = new TracePainterLine();
        ITracePainter pointPainter = new TracePainterDisc();
        /*
        IErrorBarPolicy errorBarPolicy = new ErrorBarPolicyRelative(0.0002, 0.0002);  // Percentage
        errorBarPolicy.setShowNegativeXErrors(true);
        errorBarPolicy.setShowPositiveXErrors(true);
        IErrorBarPainter errorBarPainter = new ErrorBarPainter();
        errorBarPainter.setEndPointPainter(new PointPainterDisc());
        errorBarPainter.setEndPointColor(Color.GRAY);
        errorBarPainter.setConnectionPainter(new PointPainterLine());
        errorBarPainter.setConnectionColor(Color.LIGHT_GRAY);
        errorBarPolicy.setErrorBarPainter(errorBarPainter);
        */
        Stroke thickLine = new BasicStroke(2);

        int cnt = 0;
        for(int i=0; i<nAttributes; i++)
        {
            color = StsColor.colors32[i].getColor();
            if(attributes[i] == null)
                continue;
            traces[cnt] = createDataset(i, times[i], attributes[i]);
            traces[cnt].addPropertyChangeListener(ITrace2D.PROPERTY_VISIBLE, this);
            traces[cnt].setColor(color);
            traces[cnt].setVisible(visible[i]);
            traces[cnt].setTracePainter(linePainter);
            traces[cnt].addTracePainter(pointPainter);
            traces[cnt].setStroke(thickLine);
            //traces[cnt].setErrorBarPolicy(errorBarPolicy);
            chart1.addTrace(traces[cnt]);
            cnt++;
        }
        configureAxis(chart1);
        changePolicy(chart1);
    }

    public void configureAxis(Chart2D chart1)
    {
        IAxis axis = chart1.getAxisX();
        axis.setPaintGrid(true);
        axis.setFormatter(new LabelFormatterDate(new SimpleDateFormat("kk:mm:ss")));

        axis = chart1.getAxisY();
        axis.setRange(new Range(0.0, nAttributes));
        NumberFormat numFormat = NumberFormat.getNumberInstance();
        numFormat.setMaximumFractionDigits(2);
        numFormat.setMinimumFractionDigits(0);
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
   /*
    public void setCurrentTrace()
    {
    	if(chart == null) return;
    	for(int i=0; i<chart.getTraces().size(); i++)
    	{
    		if(((Trace2DLtd)chart.getTraces().toArray()[i]).isVisible())
    			chart.setCurrentTrace(i);
    	}
    }
    */
    public void propertyChange(final PropertyChangeEvent evt)
    {
        String property = evt.getPropertyName();
        if (property.equals(ITrace2D.PROPERTY_VISIBLE))
        {
            int idx = getAttributeIndex(((ITrace2D)evt.getSource()).getLabel());
            visible[idx] = ((Boolean)evt.getNewValue()).booleanValue();
            //setCurrentTrace();
        }
        if(chart != null) chart.propertyChange(evt);
    }

    public boolean reconfigureView()
    {
    	configureAxis(chart);

    	chart.setBackground(model.getProject().getTimeSeriesBkgdColor().getColor());
    	chart.setForeground(model.getProject().getTimeSeriesFrgdColor().getColor());
    	chart.setGridColor(model.getProject().getTimeSeriesGridColor().getColor());

    	chart.getAxisX().setPaintGrid(true);
    	chart.getAxisY().setPaintGrid(false);

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
    private Trace2DLtd createDataset(int classType, long[] times, float[] attribute)
    {
        Trace2DLtd trace = null;
        trace = new Trace2DLtd(500000);
        trace.setName(names[classType]);
        for (int i = 0; i < attribute.length; i++)
        {
            trace.addPoint(times[i], attribute[i]);
        }
        return trace;
    }

    public boolean initializeMovie()
    {
        return true;
    }

    private boolean getAttributeData()
    {
        getTypeData();
        if(visible == null)
        {
            visible = new boolean[nAttributes];
            for(int i=0; i<nAttributes; i++)
                visible[i] = true;
        }
        return true;
    }

    private long[] getTimeData()
    {
        StsArrayList classList = model.getTimeDisplayableClasses();
        long[] timeArray = null;

        Iterator classIterator = classList.iterator();
        while(classIterator.hasNext())
        {
            StsClass timeClass = (StsClass)classIterator.next();
            StsObjectList objectList = timeClass.getStsObjectList();
            StsMainObject[] objects = objectList.getVisibleObjectList();
            if(objects.length != 0)
            {
                for(int i=0; i<objects.length; i++)
                {
                    StsMainTimeObject mainTime = (StsMainTimeObject)objects[i];
                    long bornDate = mainTime.getBornDateLong();
                    long deathDate = mainTime.getDeathDateLong();
                    if(bornDate > 0)
                    {
                        timeArray = StsMath.longListAddSortedValue(timeArray, bornDate);
                        if((deathDate > 0) && (deathDate > bornDate))
                            timeArray = StsMath.longListAddSortedValue(timeArray, deathDate);
                        else
                            timeArray = StsMath.longListAddSortedValue(timeArray, bornDate + 600000);
                    }
                }
                //System.out.println("Objects.length=" + objects.length);
            }
        }
        // Adding ends to the timeline so the objects at min and max are selectable
        if(timeArray != null)
        {
            timeArray = StsMath.longListAddSortedValue(timeArray, timeArray[0]-100000);
            timeArray = StsMath.longListAddSortedValue(timeArray, timeArray[timeArray.length-1]+100000);
        }
        return timeArray;
    }

    private void getTypeData()
    {
        StsArrayList classList = model.getTimeDisplayableClasses();
        Iterator classIterator = classList.iterator();
        nAttributes = classList.size();
        times = new long[nAttributes][];
        attributes = new float[nAttributes][];
        names = new String[nAttributes];
        int cnt = 0;
        while(classIterator.hasNext())
        {
            StsClass timeClass = (StsClass)classIterator.next();
            StsObjectList objectList = timeClass.getStsObjectList();
            StsObject[] objects = objectList.getElements();
            if(objects.length != 0)
            {
                int num = 0;
                for(int i=0; i<objects.length; i++)
                {
                    StsMainTimeObject mainTime = (StsMainTimeObject)objects[i];
                    long bornDate =  mainTime.getBornDateLong();
                    if(bornDate <= 0)
                        continue;

                    times[cnt] = StsMath.longListInsertValue(times[cnt], bornDate, num);
                    attributes[cnt] = StsMath.floatListInsertValue(attributes[cnt], 0.0f, num);
                    num++;
                    times[cnt] = StsMath.longListInsertValue(times[cnt], bornDate, num);
                    attributes[cnt] = StsMath.floatListInsertValue(attributes[cnt], (cnt+1), num);
                    num++;
                    if(maxTime < bornDate) maxTime = bornDate;
                    if(minTime > bornDate) minTime = bornDate;

                    long deathDate = mainTime.getDeathDateLong();
                    if((deathDate <= 0) || (deathDate < bornDate))
                        deathDate = bornDate + 600000;

                    times[cnt] = StsMath.longListInsertValue(times[cnt], deathDate, num);
                    attributes[cnt] = StsMath.floatListInsertValue(attributes[cnt], (cnt+1), num);
                    num++;
                    times[cnt] = StsMath.longListInsertValue(times[cnt], deathDate, num);
                    attributes[cnt] = StsMath.floatListInsertValue(attributes[cnt], .0f, num);
                    num++;
                    if(maxTime < deathDate) maxTime = deathDate;
                    if(minTime > deathDate) minTime = deathDate;
                }
                if(num > 0) // Found objects with a time range
                {
                    names[cnt] = timeClass.getName().substring(timeClass.getName().lastIndexOf(".Sts")+4);
                    cnt++;
                }
            }
        }
        nAttributes = cnt;
        // Add a point at both extremes to ensure proper axis sizing
        for(int i=0; i<nAttributes; i++)
        {
            if(attributes[i] == null) continue;
            times[i] = StsMath.longListInsertValue(times[i], minTime - 1, 0);
            attributes[i] = StsMath.floatListInsertValue(attributes[i], 0.0f, 0);
            times[i] = StsMath.longListInsertValue(times[i], minTime - 1, 0);
            attributes[i] = StsMath.floatListInsertValue(attributes[i], cnt, 0);
            times[i] = StsMath.longListInsertValue(times[i], maxTime + 1, attributes[i].length);
            attributes[i] = StsMath.floatListInsertValue(attributes[i], 0.0f, attributes[i].length);
            times[i] = StsMath.longListInsertValue(times[i], maxTime + 1, attributes[i].length);
            attributes[i] = StsMath.floatListInsertValue(attributes[i], cnt, attributes[i].length);
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
//        startMouseY = -1;
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
//                    endMouseY = glMouseY;
//                    endMouseY = checkYBounds(endMouseY);
//                    setYValue(endMouseY);
                    break;
                case BOTH:
                    endMouseX = mouseX;
                    endMouseX = checkXBounds(endMouseX);
                    setTime(endMouseX);
                    project.setProjectTime(endMouseX, true);

//                    endMouseY = glMouseY;
//                    endMouseY = checkYBounds(endMouseY);
//                    setYValue(endMouseY);
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
                    //chart.clearXColorRange();
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
//                    chart.clearYColorRange();
//                    startMouseY = glMouseY;
//                    startMouseY = checkYBounds(startMouseY);
//                    sensor.getSensorClass().setValueRangeAndSensor(null, 0, startMouseY, startMouseY);
//                    setYValue(startMouseY);
                    break;
                case BOTH:
                    chart.clearXColorRange();
                    startMouseX = endMouseX = mouseX;
                    startMouseX = checkXBounds(startMouseX);
                    project.clearProjectTimeDuration();
                    setByMe = false;
                    setTime(startMouseX);
                    project.setProjectTime(endMouseX, true);
//                    chart.clearYColorRange();
//                    startMouseY = glMouseY;
//                    startMouseY = checkYBounds(startMouseY);
//                    sensor.getSensorClass().setValueRangeAndSensor(null, 0, startMouseY, startMouseY);
//                    setYValue(startMouseY);
                    break;
                case NONE:
                    StsMessageFiles.infoMessage("Set ranges by clicking on X or Y axis.");
                    break;
            }
        }
    }

    private byte getActiveAxis(long x, double y)
    {
        if(y <= chart.getMinY() && x >= chart.getMinX())
            return X;
        else if(x <= chart.getMinX() && y >= chart.getMinY())
            return Y;
        else if(x < chart.getMaxX() && x > chart.getMinX() && y < chart.getMaxY() && y > chart.getMinY())
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
//                    endMouseY = glMouseY;
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
//                        endMouseY = glMouseY;
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
    		changePolicy(chart);
    		return true;
    	}
    	else if(object instanceof StsSensor)
            reconfigureView();
    	else if(object instanceof StsSensorClass)
        	reconfigureView();
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

    public void setTime(long value)
    {
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
    public void setYValue(float value)
    {
        if((value < 0) || (chart == null))
            return;
        value = checkYBounds(value);
        chart.setCrosshairYValue(value);
    }
    public void windowOpened(WindowEvent e) { }
    public void windowClosing(WindowEvent e)
    {
        project.clearProjectTimeDuration();
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
    	if(windowSize == null)
    		return;   // Rogue window is getting persisted.
    	if(chart == null) // Chart has not been restored yet.
    		return;
        setSize(windowSize);
        setLocation(windowLocation);
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
        setEnabled(true);
    }

    public Frame getFrame()
    {
        return (Frame)this.getOwner();
    }

    public Class[] getViewClasses() { return new Class[] { StsViewTimeLine.class }; }
}