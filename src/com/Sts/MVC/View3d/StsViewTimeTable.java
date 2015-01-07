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
import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.StsProject;
import com.Sts.Types.StsArrayList;
import com.Sts.UI.Beans.StsDateFieldBean;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Icons.StsIcon;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsTimeTable;
import com.Sts.UI.Toolbars.StsTimeActionToolbar;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsParameters;
import info.monitorenter.gui.chart.Chart2D;
import info.monitorenter.gui.chart.IAxis;
import info.monitorenter.gui.chart.ITrace2D;
import info.monitorenter.gui.chart.ITracePainter;
import info.monitorenter.gui.chart.controls.LayoutFactory;
import info.monitorenter.gui.chart.labelformatters.LabelFormatterDate;
import info.monitorenter.gui.chart.labelformatters.LabelFormatterNumber;
import info.monitorenter.gui.chart.traces.Trace2DLtd;
import info.monitorenter.gui.chart.traces.painters.TracePainterDisc;
import info.monitorenter.gui.chart.traces.painters.TracePainterLine;
import info.monitorenter.gui.chart.views.ChartPanel;
import info.monitorenter.util.Range;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;

public class StsViewTimeTable extends JDialog
{
    static public final String viewClassname = "TimeTable_View";
    public Dimension windowSize = null;
    public boolean[] visible = null;
    public StsTimeTable timeTable = null;

    public int timeScale = DAYS;
    final static int HOURS = 0;
    final static int DAYS = 1;
    final static int WEEKS = 2;
    final static int MONTHS = 3;
    final static int YEARS = 4;

    transient long maxTime = 0l;
    transient long minTime = StsParameters.largeLong;

    static byte STATICSENSOR = 0;
    static byte DYNAMICSENSOR = 1;
    static byte LIVEWELL = 2;
    static byte[] TIME_TYPES = new byte[] {STATICSENSOR, DYNAMICSENSOR, LIVEWELL};
    static String[] TIME_TYPE_STRINGS = new String[] {"StaticSensor", "DynamicSensor", "LiveWell"};

    transient float[][] attributes = null;
    transient String[][] names = null;
    transient long[][] times = null;
    transient int nAttributes = 0;
    transient String[] types = null;
    transient long[] time = null;
    transient boolean debug = false;

    transient StsJPanel panel = new StsJPanel();
    transient JLabel timeLabel = new JLabel("Time:");
    transient StsWin3dBase window = null;
    transient StsModel model = null;
    transient StsProject project = null;
    transient long processedTime = 0l;

    static final String[] invisibleAttributes = new String[] {"X", "Y", "Z", "DEPTH" };

    static public final String viewName = "TimeTable_View";

    /**
     * Default constructor
     */
    public StsViewTimeTable()
    {
        System.out.println("Default Time Table View Constructor.");
    }

    public StsViewTimeTable(StsWin3dBase window)
    {
        super(window, "Time Table Plot", false);
        try
        {
            this.window = window;
            this.model = window.getModel();
            this.project = model.getProject();
            buildDialog();
            setSize(800,(nAttributes+1)*16+50);
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }
    public StsViewTimeTable(StsWin3dBase window, boolean child)
    {
        super(window, "Time Table Plot", false);
        try
        {
            this.window = window;
            this.model = window.model;
            this.project = model.getProject();
            setSize(800,(nAttributes+1)*16+50);
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
        timeTable = createTable();
        buildTraces();

        panel.gbc.fill = panel.gbc.HORIZONTAL;
        panel.add(timeTable);
        processedTime = maxTime;        
    }

    void jbInit() throws Exception
    {
        add(panel);
    }

    private StsTimeTable createTable()
    {
        StsTimeTable table = new StsTimeTable(model, nAttributes,numOfColumnsInRange(), false);
        return table;
    }

    private void buildTraces()
    {
        ArrayList<Object[]> columndata = null;
        long millsPerColumn = numMillisecondsPerColumn();

        addRowNames();
        for(int i=0; i<numOfColumnsInRange(); i++)
        {
            columndata = new ArrayList<Object[]>();
            long startTime = minTime + (i*millsPerColumn);
            long endTime = startTime + millsPerColumn;
            for(int j=0; j<nAttributes; j++)
            {
                boolean found = false;
                for(int k=0; k<times[j].length; k=k+2)
                {
                    if((startTime <= times[j][k+1]) && (times[j][k] <= endTime))
                    {
                        columndata.add(new Object[]{names[j][k],StsIcon.createIcon("Timeline0" + (j+1) + ".gif"), types[j]});
                        found = true;
                        break;
                    }
                }
                if(!found)
                    columndata.add(new Object[]{"Idle "+ StsDateFieldBean.convertToString(startTime) + "-" + StsDateFieldBean.convertToString(endTime),StsIcon.createIcon("Timeline00.gif")});
            }
		    timeTable.addColumn(new Date(startTime), true, false, columndata);

            //timeTable.addColumn(null, columndata);
        }
        timeTable.setColumnWidths();
    }

    private void addRowNames()
    {
        ArrayList<Object[]> columndata = new ArrayList<Object[]>();
        columndata.add(new Object[] {"Time", null});
        for(int i=0; i<nAttributes; i++)
            columndata.add(new Object[]{types[i],null});
        timeTable.addColumn(null, columndata);
    }
    
    private int numOfColumnsInRange()
    {
        long timeRangeMs = maxTime - minTime;
        switch(timeScale)
        {
           case HOURS:
               return (int)(((timeRangeMs / 1000)/60)/60);
           case DAYS:
               return (int)((((timeRangeMs / 1000)/60)/60)/24);
           case WEEKS:
               return (int)(((((timeRangeMs / 1000)/60)/60)/24)/7);
           case MONTHS:
               return (int)(((((timeRangeMs / 1000)/60)/60)/24)/30);
           case YEARS:
               return (int)(((((timeRangeMs / 1000)/60)/60)/24)/365);
        }
        return 0;
    }

    private long numMillisecondsPerColumn()
    {
        switch(timeScale)
        {
           case HOURS:
               return 1000l*60l*60l;
           case DAYS:
               return 1000l*60l*60l*24l;
           case WEEKS:
               return 1000l*60l*60l*24l*7l;
           case MONTHS:
               return 1000l*60l*60l*24l*30l;
           case YEARS:
               return 1000l*60l*60l*24l*365l;
        }
        return 0;
    }

    public int getAttributeIndex(String name)
    {
        for(int i=0; i<types.length; i++)
        {
            if(name.equalsIgnoreCase(types[i]))
                return i;
        }
        return -1;
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

    private void getTypeData()
    {
        StsArrayList classList = model.getTimeDisplayableClasses();
        Iterator classIterator = classList.iterator();
        nAttributes = classList.size();
        times = new long[nAttributes][];
        names = new String[nAttributes][];
        attributes = new float[nAttributes][];
        types = new String[nAttributes];
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

                    long deathDate = mainTime.getDeathDateLong();
                    if((deathDate <= 0) || (deathDate < bornDate))
                        deathDate = bornDate + 600000;
                    String description = mainTime.getName() + "\n   Started:" + StsDateFieldBean.convertToString(bornDate) + "\n   End: " + StsDateFieldBean.convertToString(deathDate);

                    times[cnt] = StsMath.longListInsertValue(times[cnt], bornDate, num);
                    names[cnt] = StsMath.stringListInsertValue(names[cnt], mainTime.getName(), num);
                    attributes[cnt] = StsMath.floatListInsertValue(attributes[cnt], (cnt+1), num);
                    num++;
                    if(maxTime < bornDate) maxTime = bornDate;
                    if(minTime > bornDate) minTime = bornDate;


                    times[cnt] = StsMath.longListInsertValue(times[cnt], deathDate, num);
                    names[cnt] = StsMath.stringListInsertValue(names[cnt], mainTime.getName(), num);
                    attributes[cnt] = StsMath.floatListInsertValue(attributes[cnt], (cnt+1), num);
                    num++;
                    if(maxTime < deathDate) maxTime = deathDate;
                    if(minTime > deathDate) minTime = deathDate;
                }
                if(num > 0) // Found objects with a time range
                {
                    types[cnt] = timeClass.getName().substring(timeClass.getName().lastIndexOf(".Sts")+4);
                    cnt++;
                }
            }
        }
        nAttributes = cnt;
    }

    public boolean viewObjectChanged(Object source, Object object)
    {
        return true;
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
    }

    public void restoreWindowState()
    {
    	if(windowSize == null)
    		return;   // Rogue window is getting persisted.
        setEnabled(true);
    }

    public Frame getFrame()
    {
        return (Frame)this.getOwner();
    }

    public Class[] getViewClasses() { return new Class[] { StsViewTimeTable.class }; }
}