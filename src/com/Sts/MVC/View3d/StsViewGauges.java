package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Base class used to present two-dimensional data.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.DB.StsSerializable;
import com.Sts.DBTypes.*;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.StsProject;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.StsListSelectionDialog;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsRadioButtonDialog;
import com.Sts.UI.StsRoundGaugePanel;
import com.Sts.UI.Toolbars.StsTimeActionToolbar;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsParameters;
import com.Sts.Utilities.StsToolkit;
import info.monitorenter.gui.chart.ITrace2D;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Iterator;

public class StsViewGauges extends JDialog implements StsSerializable, WindowListener, ComponentListener
{
    transient StsWin3dBase window = null;
    transient StsModel model = null;
    transient StsProject project = null;
    transient boolean debug = false;

    public byte plotType = SENSORS;
    public Object selectedObject = null;
    public Object[] selectedCurves = null;
    public Point windowLocation = null;
    public Dimension windowSize = null;
    public int currentSensor = 0;

    transient public StsRoundGaugePanel[] gauges = null;
    static public byte SENSORS = 0;
    static public byte LIVEWELLS = 1;

    static public final String viewName = "Gauge_View";

    /**
     * Default constructor
     */
    public StsViewGauges()
    {
        System.out.println("Default Gauge View Constructor.");
    }

    public StsViewGauges(StsWin3dBase window)
    {
        super(window, "Gauge Panel", false);
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
            if(activeObjects == null) return;
            if(activeObjects.length < 1) return;

            selectedObject = StsListSelectionDialog.getSingleSelectFromListDialog(window,"Selection", null, activeObjects);
            if(selectedObject == null)
                return;

            // Determine the attributes to build gauge
            if((selectedObject instanceof StsDynamicSensor) || (selectedObject instanceof StsStaticSensor))
            {
                StsTimeCurve[] curves = ((StsSensor)selectedObject).getPropertyCurvesExcludingXYZ();
                if(curves == null)
                {
                    new StsMessage(window, StsMessage.WARNING, "No attributes available for selected sensor.");
                    return;
                }
                if(curves.length == 0)
                {
                    new StsMessage(window, StsMessage.WARNING, "No attributes available for selected sensor.");
                    return;
                }
                selectedCurves = StsListSelectionDialog.getMultiSelectFromListDialog(window,"Curve Selection", null, curves);
            }
            else
            {
                StsObjectRefList logs = ((StsLiveWell)selectedObject).getLogCurves();
                if(logs.getSize() == 0)
                {
                    new StsMessage(window, StsMessage.WARNING, "No logs available for selected well.");
                    return;
                }
                selectedCurves = StsListSelectionDialog.getMultiSelectFromListDialog(window,"Curve Selection", null, logs.getElements());
            }

            buildDialog();
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
        currentSensor = 0;
        start();
    }

    public Object[] getSelectedCurves() { return selectedCurves; }

    private StsObject[] getActiveSensorObjects(StsObject[] objects)
    {
        StsObject[] newObjects = new StsObject[objects.length];
        int cnt = 0;
        for(int i=0; i<objects.length; i++)
        {
            if(((StsSensor)objects[i]).getNumValues() > 0)
                newObjects[cnt++] = objects[i];
        }
        return (StsObject[]) StsMath.trimArray(newObjects, cnt);
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

    public void buildDialog()
    {
        try
        {
            jbInit();
            initialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
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

    }

    void jbInit() throws Exception
    {
        addWindowListener(this);
        StsGroupBox gaugeGroup = new StsGroupBox();
        int count = 0;
        gauges = new StsRoundGaugePanel[selectedCurves.length];
        for(int k=0; k<selectedCurves.length; k++)
        {
            if((selectedObject instanceof StsDynamicSensor) || (selectedObject instanceof StsStaticSensor))
            {
                StsTimeCurve curve = (StsTimeCurve)selectedCurves[k];
                int index = ((StsSensor)selectedObject).getTimeCurves().getIndex(curve);
                if(curve.getCurveMin() == curve.getCurveMax())
                {
                    gauges[k] = null;
                    continue;
                }
                gauges[k] = new StsRoundGaugePanel(curve.getName(), curve.getCurveMin(), curve.getCurveMax(), 90, StsColor.colors32[index].getColor());
            }
            else
            {
                StsLogCurve curve = (StsLogCurve)selectedCurves[k];
                int index = ((StsLiveWell)selectedObject).getLogCurves().getIndex(curve);
                if(curve.getMinValue() == curve.getMaxValue())
                {
                    gauges[k] = null;
                    continue;
                }
                gauges[k] = new StsRoundGaugePanel(curve.getName(), curve.getMinValue(), curve.getMaxValue(), 90, StsColor.colors32[index].getColor());
            }

            gaugeGroup.gbc.fill = gaugeGroup.gbc.BOTH;
            gaugeGroup.addEndRow(gauges[k]);
            count++;
        }
        this.setContentPane(gaugeGroup);
        this.addComponentListener(this);
        setSize(250, count*220);
    }

    public void clearView()
    {
    }
 
    public void resetView(StsProject project)
    {
    }

    public boolean viewObjectChanged(Object source, Object object)
    {

    	if(object instanceof StsTimeActionToolbar)
    	{
    		return true;
    	}
        // Possibly reset the min and max of the gauges.
    	else if(object instanceof StsSensor)
            ;

    	else if(object instanceof StsSensorClass)
        	;
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
                if(value < 0)
                    return;
                if(gauges == null)
                    return;
                for(int i=0; i<gauges.length; i++)
                {
                    if(gauges[i] != null)
                    {

                        if((selectedObject instanceof StsDynamicSensor) || (selectedObject instanceof StsStaticSensor))
                        {
                            StsTimeCurve curve = (StsTimeCurve)selectedCurves[i];
                            gauges[i].setValue(curve.getValueAt(value));
                        }
                        else
                        {
                            StsLogCurve curve = (StsLogCurve)selectedCurves[i];
                            StsTimeVector timeVector = ((StsLiveWell)selectedObject).getTimeVector();
                            if(timeVector == null)
                                return;
                            float indexF = timeVector.getValuesArray().length - 1;
                            if(v < timeVector.getMaxValue())
                                indexF = timeVector.getIndexF(v);
                            float depth = curve.getMDepthVector().getValue(indexF);
                            if(depth != StsParameters.largeFloat)
                            {
                                float val = curve.getInterpolatedValue(curve.getMDepthVector().getFloats(),curve.getValueVector().getFloats(), depth);
                                //System.out.println("Gauge[" + i + " ] Value=" + val);
                                gauges[i].setValue(val);
                            }
                        }
                    }
                }
            }
        };
        StsToolkit.runWaitOnEventThread(runnable);
    }

    public void windowOpened(WindowEvent e) { }
    public void windowClosing(WindowEvent e)
    {
        project.clearProjectTimeDuration();
        if(plotType == SENSORS)
            ((StsSensor)selectedObject).getSensorClass().setValueRangeAndSensor(null, 0, 0, 0);
        model.deleteAuxGaugeWindow(this);
    }
    public void windowClosed(WindowEvent e) {}
    public void windowIconified(WindowEvent e) {}
    public void windowDeiconified(WindowEvent e) {}
    public void windowActivated(WindowEvent e) {}
    public void windowDeactivated(WindowEvent e) {}
    public void componentResized(ComponentEvent e)
    {
        windowSize = new Dimension(getWidth(), getHeight());
    }
    public void componentMoved(ComponentEvent e)
    {
        windowLocation = new Point(getX(), getY());
    }
    public void componentHidden(ComponentEvent e) {}
    public void componentShown(ComponentEvent e) {}

    public void restoreWindowState()
    {
        setEnabled(true);
        if(windowSize == null)
            setSize(200,200);
        else
        {
            setSize(windowSize);
            setLocation(windowLocation);
        }
    }

    public Frame getFrame()
    {
        return (Frame)this.getOwner();
    }

    public Class[] getViewClasses() { return new Class[] { StsViewGauges.class }; }
}