package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Sounds.StsSound;
import com.Sts.Utilities.*;

import java.util.*;

public class StsWellClass extends StsLineClass implements StsSerializable, StsClassTimeDisplayable, StsClassCursorDisplayable, StsClassViewSelectable
{
//   StsWellViewsModel wellViewsModel;
    public boolean displayMarkers = true;
    public boolean displayIn2d = true;
    public boolean displayPerfMarkers = true;
    public boolean displayEquipMarkers = true;
    public boolean displayFmiMarkers = true;
    public boolean enableTime = true;
    private boolean enableSound = true;


    private boolean defaultIsDrawingCurtain = false;
    private boolean defaultIsDrawZones = false;
    private String defaultPerfSound = StsSound.IMPLOSION;
    private StsColor defaultWellColor = new StsColor(StsColor.BLUE);

    public StsLogProperties defaultLogProperties = null;
    public StsMarkerProperties defaultMarkerProperties = null;

    private String logTypeDisplay3dLeft = LOG_TYPE_NONE;
    private String logTypeDisplay3dRight = LOG_TYPE_NONE;
    /** width in pixels of log curve isVisible in 3d */
    private int logCurveDisplayWidth = 50;
    /** widht of the log curve line */
    private int LogCurveLineWidth = 1;
    /** width of line showing log curve as a colored stripe */
    private int logLineDisplayWidth = 10;
    /** indicates how logs are drawn in 3d: see displayTypeStrings */
    protected byte displayType = 3;
    private float fmiScale = 25.0f;
    private float curtainStep = 10.0f;
    private float perfScale = 25.0f;
    private float equipScale = 20.0f;
    
    private String labelFormatAsString = new String("###0.0");

    static public final String LOG_TYPE_NONE = "None";

    static final String DISPLAY_NONE = "None";
    static final String DISPLAY_CURVE = "Curve";
    static final String DISPLAY_CURTAIN = "Curtain";
    static final String DISPLAY_LINE = "Line";

    static final String[] displayTypeStrings = new String[] { DISPLAY_NONE, DISPLAY_CURVE, DISPLAY_CURTAIN, DISPLAY_LINE };

    static public final long serialVersionUID = 1L;

    public StsWellClass()
    {
        userName = "Wells";
    }

    public void initializeDisplayFields()
    {
        StsFloatFieldBean perfScaleBean = new StsFloatFieldBean(this,"perfScale", 2.0f, 50.0f, "Perforation Scale:", true);
        perfScaleBean.setRangeFixStep(2.0f, 50.0f, 1.0f);
        StsFloatFieldBean equipScaleBean = new StsFloatFieldBean(this,"equipmentScale", 2.0f, 50.0f, "Equipment Scale:", true);
        equipScaleBean.setRangeFixStep(2.0f, 50.0f, 1.0f);
        displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(this, "displayNames", "Names"),
            new StsBooleanFieldBean(this, "displayLines", "Show Path"),
            new StsBooleanFieldBean(this, "enableTime", "Enable Time"),
            new StsBooleanFieldBean(this, "enableSound", "Enable Sound"),
            new StsBooleanFieldBean(this, "displayIn2d", "Show Markers in 2D"),            
            new StsBooleanFieldBean(this, "displayMarkers", "Show Geologic Markers"),
            new StsBooleanFieldBean(this, "displayPerfMarkers", "Show Perforations"),
            new StsBooleanFieldBean(this, "displayEquipMarkers", "Show Equipment"),
            new StsBooleanFieldBean(this, "displayFmiMarkers", "Show FMI Markers"),  
            new StsFloatFieldBean(this,"fmiScale", true, "FMI Scale:"),
            perfScaleBean, equipScaleBean,
            new StsFloatFieldBean(this,"curtainStep", true, "Curtain Step:"),
            new StsStringFieldBean(this, "labelFormatAsString", "Label Format:"),
            new StsComboBoxFieldBean(this, "logTypeDisplay3dLeft", "Left Log Type:", "logTypeDisplayList"),
            new StsComboBoxFieldBean(this, "logTypeDisplay3dRight", "Right Log Type:", "logTypeDisplayList"),
            new StsComboBoxFieldBean(this, "displayTypeString", "Select log type display", displayTypeStrings),
            new StsIntFieldBean(this, "logCurveDisplayWidth", true, "Log display width, pixels"),
            new StsIntFieldBean(this, "LogCurveLineWidth", 1, 4, "Log line width, pixels", true),
            new StsIntFieldBean(this, "logLineDisplayWidth", true, "Width of colored log line, pixels")
        };
    }

    public void initializeDefaultFields()
    {
        defaultFields = new StsFieldBean[]
        {
                new StsBooleanFieldBean(this, "defaultIsDrawingCurtain", "Show Curtain"),
                new StsBooleanFieldBean(this, "defaultIsDrawZones", "Show Zones"),
                new StsComboBoxFieldBean(this, "defaultPerfSound", "Perforation Sound:", StsSound.sounds),
                new StsColorComboBoxFieldBean(this, "defaultWellColor", "Well Color:", StsColor.colors8),
        };
    }

    /** Called by logDisplayBeans to get current list of logCurveTypes to display. */
    private String[] getLogTypeDisplayList()
    {
        String[] logCurveTypeStrings = StsLogCurveType.getLogCurveTypeStrings(currentModel);
        String[] timeLogCurveTypeStrings =  StsTimeLogCurveType.getLogCurveTypeStrings();
        logCurveTypeStrings = (String[]) StsMath.arrayAddArray(logCurveTypeStrings, timeLogCurveTypeStrings);
        int nTypes = logCurveTypeStrings.length;
        String[] selectionStrings = new String[nTypes + 1];
        selectionStrings[0] = LOG_TYPE_NONE;
        System.arraycopy(logCurveTypeStrings, 0, selectionStrings, 1, nTypes);
        return selectionStrings;
    }

    public String getLogTypeDisplay3dLeft() { return logTypeDisplay3dLeft; }

    public void setLogTypeDisplay3dLeft(String type)
    {
        if(this.logTypeDisplay3dLeft == type) return;
        this.logTypeDisplay3dLeft = type;
        StsObject[] wellObjects = getElements();
        for(int n = 0; n < wellObjects.length; n++)
            ((StsWell)wellObjects[n]).setLogTypeDisplay3dLeft(type);
        currentModel.repaintViews(StsView3d.class);
    }

    public String getLogTypeDisplay3dRight() { return logTypeDisplay3dRight; }

    public void setLogTypeDisplay3dRight(String type)
    {
        if(this.logTypeDisplay3dRight == type) return;
        this.logTypeDisplay3dRight = type;
        StsObject[] wellObjects = getElements();
        for(int n = 0; n < wellObjects.length; n++)
            ((StsWell)wellObjects[n]).setLogTypeDisplay3dRight(type);
        currentModel.repaintViews(StsView3d.class);
    }

    public int getLogCurveDisplayWidth() { return logCurveDisplayWidth; }

    public void setLogCurveDisplayWidth(int width)
    {
        if(logCurveDisplayWidth == width) return;
        logCurveDisplayWidth = width;
        currentModel.repaintViews(StsView3d.class);
    }

    public int getLogLineDisplayWidth() { return logLineDisplayWidth; }

    public void setLogLineDisplayWidth(int width)
    {
        if(logLineDisplayWidth == width) return;
        logLineDisplayWidth = width;
        currentModel.repaintViews(StsView3d.class);
    }

    public void setDisplayTypeString(String typeString)
    {
        for(int n = 0; n < displayTypeStrings.length; n++)
            if(displayTypeStrings[n] == typeString)
                setDisplayType((byte)n);
    }

    private void setDisplayType(byte type)
    {
        if(displayType == type) return;
        displayType = type;
        currentModel.viewObjectRepaint(this, this);
    }

    public String getDisplayTypeString() { return displayTypeStrings[displayType]; }
    public float getEquipmentScale() { return equipScale; }
    public void setEquipmentScale(float scale)
    {
        if(equipScale == scale) return;
        equipScale = scale;
//        setDisplayField("equipScale", equipScale);
        currentModel.repaintViews(StsView3d.class);
    }

    public float getPerfScale() { return perfScale; }
    public void setPerfScale(float scale)
    {
        if(perfScale == scale) return;
        perfScale = scale;
 //      setDisplayField("perfScale", perfScale);        
        currentModel.repaintViews(StsView3d.class);
    }
    public float getCurtainStep() { return curtainStep; }
    public void setCurtainStep(float step)
    {
        if(curtainStep == step) return;
        curtainStep = step;
 //      setDisplayField("curtainStep", curtainStep);
        
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsWell well = (StsWell)iter.next();
            well.reconfigureCurtainOffsetBeans();
        }        
        currentModel.repaintViews(StsView3d.class);
    }
    public float getFmiScale() { return fmiScale; }
    public void setFmiScale(float scale)
    {
        if(fmiScale == scale) return;
        fmiScale = scale;
 //      setDisplayField("fmiScale", fmiScale);
        currentModel.repaintViews(StsView3d.class);
    }
    
    public int getLogCurveLineWidth() { return LogCurveLineWidth; }

    public void setLogCurveLineWidth(int width)
    {
        if(LogCurveLineWidth == width) return;
        LogCurveLineWidth = width;
        currentModel.repaintViews(StsView3d.class);
    }

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsWell well = (StsWell)iter.next();
            if(displayLines)
                well.display(glPanel3d, displayNames, displayMarkers, displayPerfMarkers, displayFmiMarkers, this);
            if(well.wellViewModel == null) continue;
            if(well.wellViewModel.isVisible)
                well.wellViewModel.display(glPanel3d);
        }
    }

    public void displayTimeClass(StsGLPanel3d glPanel3d, long time)
    {
        if(!enableTime)
        {
            displayClass(glPanel3d);
            return;
        }

        // Possibly add list here with time settings for all objects in class to sort and determine which to display
        // currently the instance will decide whether to display itself.

        // Cycle through all well instances.
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsWell well = (StsWell)iter.next();
            if(!well.isAlive(time))
                continue;
            if(displayLines)
                well.display(glPanel3d, displayNames, displayMarkers, displayPerfMarkers, displayFmiMarkers, this);
            if(well.wellViewModel == null) continue;
            if(well.wellViewModel.isVisible)
                well.wellViewModel.display(glPanel3d);
        }
    }
         
    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging)
    {
    }

    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsWell well = (StsWell)iter.next();
            if(displayLines)
            {
                if(!currentModel.getProject().isDepth())
                {
                    if(well.isInTime())
                        well.display2d(glPanel3d, displayNames, dirNo, dirCoordinate, axesFlipped,
                                xAxisReversed, yAxisReversed, displayIn2d,
                                displayMarkers, displayPerfMarkers, displayFmiMarkers);
                }
                else
                    well.display2d(glPanel3d, displayNames, dirNo, dirCoordinate, axesFlipped,
                            xAxisReversed, yAxisReversed, displayIn2d,
                            displayMarkers, displayPerfMarkers, displayFmiMarkers);
            }
        }
    }

    public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain)
    {
    }

    public void displayLogMarkerProperties()
    {
        StsDialogFace[] dialogs = {defaultLogProperties, defaultMarkerProperties};
        new StsOkApplyCancelDialog(currentModel.win3d, dialogs , "Edit Log & Marker Properties", false);
    }

    public boolean projectInitialize(StsModel model)
    {
        defaultLogProperties = new StsLogProperties(model, this, "defaultLogProperties");
        defaultMarkerProperties = new StsMarkerProperties(model, this, "defaultMarkerProperties");
        return true;
    }

    public void setDisplayNames(boolean displayNames)
    {
        if(this.displayNames == displayNames) return;
        this.displayNames = displayNames;
        /*
        for(int n = 0; n < getSize(); n++)
        {
            StsWell well = (StsWell)getElement(n);
            well.setDrawLabels(displayNames);
        } 
		*/      
 //      setDisplayField("displayNames", displayNames);
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplayNames() 
    { 
    	return displayNames; 
    }

    public void setEnableTime(boolean enable)
    {
        if(this.enableTime == enable) return;
        this.enableTime = enable;
 //      setDisplayField("enableTime", enableTime);
        currentModel.win3dDisplayAll();
    }

    public boolean getEnableTime() { return enableTime; }

    public void setDisplayLines(boolean displayLines)
    {
        if(this.displayLines == displayLines) return;
        this.displayLines = displayLines;
 //      setDisplayField("displayLines", displayLines);
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplayLines() { return displayLines; }
    
    public void setDisplayEquipMarkers(boolean displayEquipMarkers)
    {
        if(this.displayEquipMarkers == displayEquipMarkers) return;
        this.displayEquipMarkers = displayEquipMarkers;
        currentModel.win3dDisplayAll();
    }
    public boolean getDisplayEquipMarkers() { return displayEquipMarkers; }
    public void setDisplayPerfMarkers(boolean displayPerfMarkers)
    {
        if(this.displayPerfMarkers == displayPerfMarkers) return;
        this.displayPerfMarkers = displayPerfMarkers;
 //      setDisplayField("displayPerfMarkers", displayPerfMarkers);
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplayPerfMarkers() { return displayPerfMarkers; }
    
    public void setDisplayFmiMarkers(boolean displayFmiMarkers)
    {
        if(this.displayFmiMarkers == displayFmiMarkers) return;
        this.displayFmiMarkers = displayFmiMarkers;
 //      setDisplayField("displayFmiMarkers", displayFmiMarkers);
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplayFmiMarkers() { return displayFmiMarkers; }
        
    public void setDisplayMarkers(boolean displayMarkers)
    {
        if(this.displayMarkers == displayMarkers) return;
        this.displayMarkers = displayMarkers;
 //      setDisplayField("displayMarkers", displayMarkers);
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplayIn2d() { return displayIn2d; }
    
    public void setDisplayIn2d(boolean display2d)
    {
        if(this.displayIn2d == display2d) return;
        this.displayIn2d = display2d;
 //      setDisplayField("displayIn2d", displayIn2d);
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplayMarkers() { return displayMarkers; }

    public void setDefaultIsDrawZones(boolean defaultIsDrawZones)
    {
        if(this.defaultIsDrawZones == defaultIsDrawZones) return;
        this.defaultIsDrawZones = defaultIsDrawZones;
    }

    public boolean getDefaultIsDrawZones() { return defaultIsDrawZones; }

    public void setDefaultIsDrawingCurtain(boolean defaultIsDrawingCurtain)
    {
        if(this.defaultIsDrawingCurtain == defaultIsDrawingCurtain) return;
        this.defaultIsDrawingCurtain = defaultIsDrawingCurtain;
    }

    public boolean getDefaultIsDrawingCurtain() { return defaultIsDrawingCurtain; }

    public void setLabelFormatAsString(String labelFormat)
    {
        if(this.labelFormatAsString == labelFormat) return;
        this.labelFormatAsString = labelFormat;
//        setDisplayField("labelFormatAsString", labelFormatAsString);
        currentModel.win3dDisplayAll();
    }

    public String getLabelFormatAsString() { return labelFormatAsString; }

    //TODO should change ColorListBean to use StsColors rather than Java colors
    public void setDefaultWellColor(StsColor color)
    {
        if(defaultWellColor.equals(color)) return;
        defaultWellColor = color;
    }
    //TODO should change ColorListBean to use StsColors rather than Java colors 
    public StsColor getDefaultWellColor() { return defaultWellColor; }

    public StsColor getDefaultWellStsColor() { return defaultWellColor; }

    public boolean viewObjectChanged(Object source, Object object)
    {
        for(int n = 0; n < getSize(); n++)
        {
            StsWellFrameViewModel wellViewModel = ((StsWell)getElement(n)).wellViewModel;
            if(wellViewModel == null) continue;
            if(wellViewModel.isVisible)
                wellViewModel.viewObjectChanged(source, object);
        }

        return false;
    }

    public boolean viewObjectRepaint(Object source, Object object)
    {
        for(int n = 0; n < getSize(); n++)
        {
            StsWellFrameViewModel wellViewModel = ((StsWell)getElement(n)).wellViewModel;
            if(wellViewModel == null) continue;
            if(wellViewModel.isVisible)
                wellViewModel.viewObjectRepaint(source, object);
        }

        return false;
    }

    public void closeWindows()
    {
        for(int n = 0; n < getSize(); n++)
        {
            StsWellFrameViewModel wellViewModel = ((StsWell)getElement(n)).wellViewModel;
            if(wellViewModel == null) continue;
            if(wellViewModel.isVisible)
                wellViewModel.closeWindow();
        }
    }

    public void initializePropertyFields()
    {
        propertyFields = new StsFieldBean[]
        {
                new StsButtonFieldBean("Log & Marker Properties", "Edit Log & Marker properties.", this, "displayLogMarkerProperties")
        };
    }

    public StsColor getDefaultColor(int index)
    {
        if(defaultMarkerProperties != null)
            return defaultMarkerProperties.getColor(index);
        else
            return StsColor.colors32[index];
    }
    public void setEnableSound(boolean enable)
    {
        if(this.enableSound == enable) return;
        this.enableSound = enable;
        currentModel.win3dDisplayAll();
    }

    public boolean getEnableSound() { return enableSound; }
    public void setDefaultPerfSound(String sound)
    {
        if (this.defaultPerfSound.equalsIgnoreCase(sound)) return;
        StsSound.play(sound);
        this.defaultPerfSound = sound;
    }


    public String getDefaultPerfSound() { return defaultPerfSound; }
/*
    public StsWell selectWell() throws StsException
    {
        // retrieve a well
        StsSelectStsObjects selector = StsSelectStsObjects.constructor(this,
            StsWell.class, name, "Select a well:", true);
        if (selector == null) return null;
        StsWell well = (StsWell) selector.selectObject();
        if (well == null)
        {
            logMessage("No well selected.");
            return null;
        }
        return well;
    }
*/
}