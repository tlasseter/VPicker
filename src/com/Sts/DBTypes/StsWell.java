//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Actions.Import.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.StsComboBoxToolbar;
import com.Sts.Utilities.*;
import com.Sts.WorkflowPlugIn.StsWorkflowPlugIn;
import com.magician.fonts.*;

import javax.media.opengl.*;
import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.text.*;
import java.util.*;

public class StsWell extends StsLine implements StsTimeEnabled, StsSelectable, StsTreeObjectI, StsMonitorable, StsViewSelectable
{
    // instance fields
    protected StsObjectRefList logCurves = null;
    protected StsObjectRefList lithZones = null;
    protected StsObjectRefList markers = null;
    protected StsLogCurve[] leftDisplayLogCurves;
    protected StsLogCurve[] rightDisplayLogCurves;
    /** StsSurfaceVertex at intersection of well, surface, block */
    //    protected StsObjectRefList surfaceVertices = null;
    protected boolean drawMarkers = true;
    protected boolean drawSurfaceMarkers = true;
    protected boolean drawPerfMarkers = true;
    protected boolean drawFMIMarkers = true;
    protected boolean drawEquipmentMarkers = true;
    protected boolean hasVsp = false;  // used in vsp display on well window; should be String vspName

    transient public StsWellFrameViewModel wellViewModel;

    transient protected StsSeismicCurtain seismicCurtain = null;

    protected String drawLabelString = NO_LABEL;
    protected float labelInterval = 100.0f;
    /** Alarms */
    protected StsAlarm[] alarms = null;

    static final String WELL_DEV_PREFIX = StsLogVector.WELL_DEV_PREFIX;
    static final String WELL_LOG_PREFIX = StsLogVector.WELL_LOG_PREFIX;
    static final String WELL_TD_PREFIX = StsLogVector.WELL_TD_PREFIX;
    static final String WELL_REF_PREFIX = "well-ref";

    static final int WEBJARFILE = 0;
    static final int JARFILE = 1;
    static final int BINARYFILES = 2;
    static final int ASCIIFILES = 3;
    static final int LASFILES = 4;
    static final int UTFILES = 5;

    static final String[] TD_BOTH_STRINGS = StsParameters.TD_BOTH_STRINGS;

    // Maintain the native units
    protected byte nativeHorizontalUnits = StsParameters.DIST_FEET;
    protected byte nativeVerticalUnits = StsParameters.DIST_FEET;

    protected float curtainOffset = 0.0f;

    // Stripped from imported well file headers
    protected String operator = "unknown";
    protected String company = "unknown";
    protected String field = "unknown";
    protected String area = "unknown";
    protected String state = "unknown";
    protected String county = "unknown";
    protected String wellNumber = "unknown";
    protected String wellLabel = "unknown";
    protected String api = "00000000000000";
    protected String uwi = "00000000000000";
    protected String date = "unknown";
    protected float kbElev = 0.0f;
    protected float elev = 0.0f;
    protected String elevDatum = "unknown";
    protected long spudDate = 0L;
    protected long completionDate = 0L;
    protected long permitDate = 0L;


    public boolean isDrawingCurtain = false;
    public boolean drawCurtainTransparent = false;

    public transient boolean hasMDepths = false;

    //	private transient float[] verticesDepths = null;
    //	private transient float[] verticesMDepths = null;
    //	private transient float[] verticesTimes = null;
    //	private transient StsPoint[] timePoints = null;

    //    private transient String exportName = null;
    //    private transient boolean exportLogData;
    //    private transient boolean exportDeviationData;
    //    private transient boolean exportSeisAttData;

    static final float maxMDepthError = 10.0f;

    static protected StsWell currentWell = null;
    static protected StsObjectPanel objectPanel = null;

    // Convenience flag copies
    static public final int NONE = StsParameters.NONE;
    static public final int STRAT = StsParameters.STRAT;
    static public final int LITH = StsParameters.LITH;

    static public final float nullValue = StsParameters.nullValue;
    static public final String TIME = StsWellKeywordIO.TIME;

    public transient DecimalFormat labelFormat = new DecimalFormat("###0.0");

    // display fields

    static public final String NO_LABEL = "No Label";
    static public final String MDEPTH_LABEL = "Measured Depth";
    static public final String DEPTH_LABEL = "TVD";
    static public final String TIME_LABEL = "Time";

    static public final String[] LABEL_STRINGS = new String[]{NO_LABEL, MDEPTH_LABEL, DEPTH_LABEL, TIME_LABEL};

    static StsDateFieldBean bornField = null;
    static StsDateFieldBean deathField = null;
    static StsFloatFieldBean curtainOffsetBean = null;
    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;

    /** default constructor */
    public StsWell() // throws StsException
    {
    }

    public StsWell(boolean persistent) // throws StsException
    {
        super(persistent);
    }

    public StsWell(String name, boolean persistent)
    {
        this(name, persistent, getWellClass().getDefaultWellStsColor());
    }

    public StsWell(String name, boolean persistent, StsColor color)
    {
        super(name, persistent);
        setZDomainOriginal(StsParameters.TD_DEPTH);
        isDrawingCurtain = getWellClass().getDefaultIsDrawingCurtain();
        drawZones = getWellClass().getDefaultIsDrawZones();
        setVerticesRotated(false);
        setStsColor(color);
    }

    static public StsWell nullWellConstructor(String name)
    {
        return new StsWell(name, false);
    }

    public void initializeColor()
    {
        if (currentModel != null)
            stsColor = new StsColor(currentModel.getProject().getDefaultWellStsColor());
        else
            stsColor = new StsColor(Color.RED);
    }

    /*
       static public FilenameSet createFilenameSet(String wellname, int type)
       {
           return new FilenameSet(wellname, type);
       }

       static public class FilenameSet
       {
           public String wellname;
           public int welltype = -1;
           public int nLogFilenames = 0;
           public String[] devFilenames;
           public String[] logFilenames = null;
           public int nRefFilenames = 0;
           public String[] refFilenames = null;
           int increment = 10;

           public boolean hasDepth = false;

           public FilenameSet(String wellname, int welltype)
           {
               this.wellname = wellname;
               this.welltype = welltype;
               devFilenames = new String[4];
               logFilenames = new String[4];
           }

           // add entries according to type; see StsLogVector for type definitions
           // types: 0 = X, 1 = Y, 2 = DEPTH, 3 = MDEPTH
           public void addBinaryFilename(String filename, String prefix, String curvename)
           {
               byte type = StsLogVector.getTypeFromString(StsWellKeywordIO.getCurveName());

               if (prefix.equals(WELL_DEV_PREFIX))
               {
                   devFilenames[type] = filename;
               }
               else if (prefix.equals(WELL_LOG_PREFIX))
               {
                   logFilenames = (String[]) StsMath.arrayAddElement(logFilenames, filename, nLogFilenames, increment);
                   nLogFilenames++;
                   if (type == StsLogVector.DEPTH || type == StsLogVector.MDEPTH)
                   {
                       hasDepth = true;
                   }
               }
           }

           public void addAsciiFilename(String filename, String prefix, int type)
           {
               if (prefix.equals(WELL_DEV_PREFIX))
               {
                   devFilenames[0] = filename;
                   welltype = type;
               }
               else if (prefix.equals(WELL_LOG_PREFIX) || prefix.equals(WELL_TD_PREFIX))
               {
                   logFilenames = (String[]) StsMath.arrayAddElement(logFilenames, filename, nLogFilenames, increment);
                   nLogFilenames++;
               }
               else if(prefix.equals(WELL_REF_PREFIX))
               {
                   refFilenames = (String[]) StsMath.arrayAddElement(refFilenames, filename, nRefFilenames, increment);
                   nRefFilenames++;
               }
           }

           public boolean binaryAlreadyInList(String filename, String prefix, String curvename)
           {
               byte type = StsLogVector.getTypeFromString(StsWellKeywordIO.getCurveName());
               if (prefix.equals(WELL_DEV_PREFIX))
               {
                   if (devFilenames[type] != null)
                   {
                       return true;
                   }
               }
               else if (prefix.equals(WELL_LOG_PREFIX))
               {
                   for (int i = 0; i < nLogFilenames; i++)
                   {
                       if (logFilenames[i].equals(filename))
                       {
                           return true;
                       }
                   }
               }
               return false;
           }

           public boolean asciiAlreadyInList(String filename, String prefix)
           {
               if (prefix.equals(WELL_DEV_PREFIX))
               {
                   if (devFilenames[0] != null)
                   {
                       return true;
                   }
                   else if (prefix.equals(WELL_LOG_PREFIX))
                   {
                       if (logFilenames[0].equals(filename))
                       {
                           return true;
                       }
                   }
               }
               return false;
           }

           public boolean isOK()
           {
               if ((welltype == ASCIIFILES) || (welltype == LASFILES) || (welltype == UTFILES))
               {
                   if (devFilenames[0] == null)
                   {
                       return false;
                   }
                   return true;
               }
               else
               {
                   // check deviation vectors
                   if (devFilenames[0] == null || devFilenames[1] == null || (devFilenames[2] == null && devFilenames[3] == null))
                   {
                       return false;
                   }

                   if (!hasDepth)
                   {
                       logFilenames = new String[0];
                       nLogFilenames = 0;
                       return false;
                   }
                   trimLogFilenames();
                   StsMath.qsort(logFilenames);
                   return true;
               }
           }

           public void trimLogFilenames()
           {
               logFilenames = (String[]) StsMath.trimArray(logFilenames, nLogFilenames);
           }
       }
    */
    public boolean isTimeEnabled()
    {
        return getWellClass().getEnableTime();
    }
    public float getTopZ()
    {
        if(lineVertices == null)
            return 0.0f;
        return getTopPoint().getZ();
    }

    public float getBotZ()
    {
        if(lineVertices == null)
           return 0.0f;
        return getBotPoint().getZ();
    }

    public double getXOrigin()
    {
        return xOrigin;
    }

    public double getYOrigin()
    {
        return yOrigin;
    }

    public float getCurtainOffset() { return curtainOffset; }

    public void setCurtainOffset(float val)
    {
        if (isVertical && val > 0.0f)
        {
            new StsMessage(currentModel.win3d, StsMessage.INFO, "Can't offset vertical well: resetting curtain offset to zero");
//            StsMessageFiles.errorMessage("Can't offset vertical well: resetting curtain offset to zero");
            curtainOffset = 0.0f;
            curtainOffsetBean.setValue(curtainOffset);
            dbFieldChanged("curtainOffset", curtainOffset);
            return;
        }

        if (curtainOffset == val) return;

        curtainOffset = val;
        dbFieldChanged("curtainOffset", curtainOffset);
        deleteSeismicCurtain();
        createCurtain();
        currentModel.viewObjectChanged(this, this);
    }

    public String getOperator() { return operator; }

    public String getField() { return field; }

    public String getArea() { return area; }

    public String getState() { return state; }

    public String getCounty() { return county; }

    public String getWellNumber() { return wellNumber; }

    public String getWellLabel() { return wellLabel; }

    public String getApi() { return api; }

    public String getUwi() { return uwi; }

    public String getDate() { return date; }

    public String getCompany() { return company; }

    public float getKbElev() { return kbElev; }

    public float getElev() { return elev; }

    public String getElevDatum() { return elevDatum; }

    public long getSpudDate() { return spudDate; }

    public String getSpudDateString()
    {
        if (spudDate == 0l)
            return "Undefined";
        else
            return currentModel.getProject().getDateFormat().format(new Date(spudDate));
    }

    public long getCompletionDate() { return completionDate; }

    public String getCompletionDateString()
    {
        if (completionDate == 0l)
            return "Undefined";
        else
            return currentModel.getProject().getDateFormat().format(new Date(completionDate));
    }

    public long getPermitDate() { return permitDate; }

    public String getPermitDateString()
    {
        if (permitDate == 0l)
            return "Undefined";
        else
            return currentModel.getProject().getDateFormat().format(new Date(permitDate));
    }

    public void setKbElev(float kbEl) { kbElev = kbEl; }

    public void setElev(float el) { elev = el; }

    public void setElevDatum(String ed) { elevDatum = ed; }

    public void setOperator(String op) { operator = op; }

    public void setField(String fld) { field = fld; }

    public void setArea(String a) { area = a; }

    public void setState(String st) { state = st; }

    public void setCounty(String cnty) { county = cnty; }

    public void setWellNumber(String num) { wellNumber = num; }

    public void setWellLabel(String label) { wellLabel = label; }

    public void setApi(String api) { this.api = api; }

    public void setUwi(String uwi) { this.uwi = uwi; }

    public void setDate(String date) { this.date = date; }

    public void setCompany(String company) { this.company = company; }

    public void setSpudDate(long spud) { this.spudDate = spud; }

    public void setCompletionDate(long complete) { this.completionDate = complete; }

    public void setPermitDate(long permit) { this.permitDate = permit; }

    public void setSpudDateString(String spud) { }

    public void setCompletionDateString(String complete) { }

    public String getNativeVerticalUnitsString()
    {
        return StsParameters.DIST_STRINGS[nativeVerticalUnits];
    }

    public byte getNativeVerticalUnits()
    {
        return nativeVerticalUnits;
    }

    public String getNativeHorizontalUnitsString()
    {
        return StsParameters.DIST_STRINGS[nativeHorizontalUnits];
    }

    public byte getNativeHorizontalUnits()
    {
        return nativeHorizontalUnits;
    }

    static public void initColors()
    {
        StsColorComboBoxFieldBean colorComboBox;
        StsSpectrum spectrum = currentModel.getSpectrum("Basic");
        colorComboBox = (StsColorComboBoxFieldBean) StsFieldBean.getBeanWithFieldName(displayFields, "stsColor");
        colorComboBox.setListItems(spectrum);
    }

    public StsFieldBean[] getDisplayFields()
    {
        if (displayFields == null)
        {
            bornField = new StsDateFieldBean(StsWell.class, "bornDate", "Born Date:");
            deathField = new StsDateFieldBean(StsWell.class, "deathDate", "Death Date:");
            curtainOffsetBean = new StsFloatFieldBean(StsWell.class, "curtainOffset", true, "CurtainOffset", true);
            curtainOffsetBean.setRangeFixStep(-1000.0, 1000.0, getWellClass().getCurtainStep());
            displayFields = new StsFieldBean[]
                {
                    new StsBooleanFieldBean(StsWell.class, "isVisible", "Enable"),
                    bornField,
                    deathField,
                    new StsBooleanFieldBean(StsWell.class, "drawZones", "Show Zones"),
                    new StsBooleanFieldBean(StsWell.class, "drawMarkers", "Show Markers"),
                    new StsBooleanFieldBean(StsWell.class, "drawPerfMarkers", "Show Perforations"),
                    new StsBooleanFieldBean(StsWell.class, "drawEquipmentMarkers", "Show Equipment"),
                    new StsBooleanFieldBean(StsWell.class, "drawCurtain", "Show Curtain"),
                    new StsBooleanFieldBean(StsWell.class, "drawCurtainTransparent", "Make Curtain Transparent"),
                    curtainOffsetBean,
                    new StsButtonFieldBean("Select left logs", "Select logs to display on left of well from dialog.", this, "logDisplay3dLeft"),
                    new StsButtonFieldBean("Select right logs", "Select logs to display on right of well from dialog.", this, "logDisplay3dRight"),
                    new StsBooleanFieldBean(StsWell.class, "highlighted", "Highlight Well"),
                    new StsComboBoxFieldBean(StsWell.class, "drawLabelString", "Label Type:", LABEL_STRINGS),
                    new StsFloatFieldBean(StsWell.class, "labelInterval", true, "Label Interval:"),
                    new StsColorComboBoxFieldBean(StsWell.class, "stsColor", "Color:", currentModel.getSpectrum("Basic").getStsColors())
                };
        }
        return displayFields;
    }

    public void reconfigureCurtainOffsetBeans()
    {
        if (curtainOffsetBean != null)
            curtainOffsetBean.setRangeFixStep(-1000.0, 1000.0, getWellClass().getCurtainStep());
    }

    public StsFieldBean[] getPropertyFields()
    {
        if(propertyFields == null)
        {
            propertyFields = new StsFieldBean[]
            {
                    new StsFloatFieldBean(StsWell.class, "topZ", false, "Min Depth:"),
                    new StsFloatFieldBean(StsWell.class, "botZ", false, "Max Depth:"),
                    new StsDoubleFieldBean(StsWell.class, "xOrigin", false, "X Origin:"),
                    new StsDoubleFieldBean(StsWell.class, "yOrigin", false, "Y Origin:"),
                    new StsStringFieldBean(StsWell.class, "operator", false, "Operator:"),
                    new StsStringFieldBean(StsWell.class, "field", false, "Field:"),
                    new StsStringFieldBean(StsWell.class, "area", false, "Area:"),
                    new StsStringFieldBean(StsWell.class, "state", false, "State:"),
                    new StsStringFieldBean(StsWell.class, "county", false, "County:"),
                    new StsStringFieldBean(StsWell.class, "wellNumber", false, "Well Number:"),
                    new StsStringFieldBean(StsWell.class, "wellLabel", false, "Well Label:"),
                    new StsStringFieldBean(StsWell.class, "api", false, "API Number:"),
                    new StsStringFieldBean(StsWell.class, "uwi", false, "UWI Number:"),
                    new StsStringFieldBean(StsWell.class, "date", false, "Date:"),
                    new StsStringFieldBean(StsWell.class, "company", false, "Company:"),
                    new StsFloatFieldBean(StsWell.class, "kbElev", false, "Kelly Elevation:"),
                    new StsFloatFieldBean(StsWell.class, "elev", false, "Elevation:"),
                    new StsStringFieldBean(StsWell.class, "elevDatum", false, "Elevation Datum:"),
                    new StsStringFieldBean(StsWell.class, "nativeHorizontalUnitsString", false, "Native Horizontal Units:"),
                    new StsStringFieldBean(StsWell.class, "nativeVerticalUnitsString", false, "Native Vertical Units:"),
                    new StsStringFieldBean(StsWell.class, "spudDateString", false, "Spud Date:"),
                    new StsStringFieldBean(StsWell.class, "completionDateString", false, "Completion Date:"),
                    new StsStringFieldBean(StsWell.class, "permitDateString", false, "Permit Date:")
            };
        }
        return propertyFields;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        getWellClass().selected(this);
    }

    static public StsWellClass getWellClass()
    {
        return (StsWellClass) currentModel.getCreateStsClass(StsWell.class);
    }

    /** for wells, measuredDepth should be read or computed when well is loaded */
    protected void addMeasuredDepth(StsPoint[] projectPoints)
    {
    }

    public boolean projectRotationAngleChanged()
    {
        if (!computePoints()) return false;
        if (markers == null) return false;
        markers.forEach("resetLocation");
        return true;
    }

    public void display2d(StsGLPanel3d glPanel3d, boolean displayName, int dirNo,
                          float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed,
                          boolean displayIn2d, boolean displayAllMarkers,
                          boolean displayPerfMarkers, boolean displayFmiMarkers)
    {
        if (!currentModel.getProject().canDisplayZDomain(zDomainSupported))
        {
            return;
        }
        if (glPanel3d == null)
        {
            return;
        }
        super.display2d(glPanel3d, displayName, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);

        if ((markers != null) && (displayIn2d))
        {
            int nMarkers = markers.getSize();
            int nPerfs = 0;
            for (int n = 0; n < nMarkers; n++)
            {
                StsWellMarker marker = (StsWellMarker) markers.getElement(n);
                boolean drawDifferent = marker.getMarker().getModelSurface() != null;
                if (marker instanceof StsPerforationMarker)
                {
                    if ((displayPerfMarkers) && (drawPerfMarkers))
                        ((StsPerforationMarker) marker).display2d(glPanel3d, dirNo, displayName, StsWell.getWellClass().getDefaultColor(nPerfs) , 1.0f);
                    nPerfs++;
                }
                else if (marker instanceof StsEquipmentMarker)
                {
                    if (drawEquipmentMarkers)
                        ((StsEquipmentMarker) marker).display2d(glPanel3d, dirNo, displayName, drawDifferent);
                }
                else if (marker instanceof StsFMIMarker)
                {
                    if (displayFmiMarkers)
                        ((StsFMIMarker) marker).display2d(glPanel3d, dirNo, displayName);
                }
                else
                {
                    if ((displayAllMarkers) && (drawMarkers))
                        marker.display2d(glPanel3d, dirNo, displayName, drawDifferent);
                }
            }
        }
    }

    public void display(StsGLPanel3d glPanel3d, boolean displayName, boolean displayAllMarkers, StsWellClass wellClass)
    {
        display(glPanel3d, displayName, displayAllMarkers, true, true, wellClass);
    }

    public void display(StsGLPanel3d glPanel3d, boolean displayName, boolean displayAllMarkers,
                        boolean displayPerfMarkers, boolean displayFmiMarkers, StsWellClass wellClass)
    {
        if (!currentModel.getProject().canDisplayZDomain(zDomainSupported)) return;
        if (glPanel3d == null) return;
        if (rotatedPoints == null) return;

        labelFormat = new DecimalFormat(getWellClass().getLabelFormatAsString());
        /*
                if(debugMarkersPrint)
                {
                    debugCheckMarkers();
                    debugMarkersPrint = false;
                }
        */
        if (isDrawingCurtain)
        {
            glPanel3d.setViewShift(glPanel3d.getGL(), 2.0f);
            if (displayName)
                super.display(glPanel3d, true, getName(), rotatedPoints);
            else
                super.display(glPanel3d, true, null, rotatedPoints);
        }
        else
        {
            if (displayName)
                super.display(glPanel3d, getName(), rotatedPoints);
            else
                super.display(glPanel3d, null, rotatedPoints);
        }
        if (isDrawingCurtain)
        {
            glPanel3d.resetViewShift(glPanel3d.getGL());
        }
        GL gl = glPanel3d.getGL();

        if (markers != null)
        {
            int nMarkers = markers.getSize();
            int nPerfs = 0;
            for (int n = 0; n < nMarkers; n++)
            {
                StsWellMarker marker = (StsWellMarker) markers.getElement(n);
                boolean drawDifferent = marker.getMarker().getModelSurface() != null;
                if (marker instanceof StsPerforationMarker)
                {
                    if ((displayPerfMarkers) && (drawPerfMarkers))
                        ((StsPerforationMarker)marker).display(glPanel3d, displayName, isDrawingCurtain, StsWell.getWellClass().getDefaultColor(nPerfs));
                    nPerfs++;
                }
                else if (marker instanceof StsEquipmentMarker)
                {
                    if (drawEquipmentMarkers)
                        marker.display(glPanel3d, displayName, isDrawingCurtain, drawDifferent);
                }
                else if (marker instanceof StsFMIMarker)
                {
                    if (displayFmiMarkers)
                        marker.display(glPanel3d, displayName, isDrawingCurtain, drawDifferent);
                }
                else
                {
                    if ((displayAllMarkers) && (drawMarkers))
                        marker.display(glPanel3d, displayName, isDrawingCurtain, drawDifferent);
                }
            }
        }

        logDisplay3d(leftDisplayLogCurves, -1f, glPanel3d);
        logDisplay3d(rightDisplayLogCurves, 0f, glPanel3d);

        if (!drawLabelString.equalsIgnoreCase(NO_LABEL) && (labelInterval >= 0.0f))
        {
            StsPoint point = null;
            float md = 0.0f;
            String label = null;
            int nLabels = (int) (getMaxMDepth() / labelInterval);

            if (isDrawingCurtain)
            {
                StsColor.BLACK.setGLColor(gl);
                glPanel3d.setViewShift(gl, 10.0f);
            }
            else
            {
                stsColor.setGLColor(gl);
                glPanel3d.setViewShift(gl, 1.0f);
            }

            GLBitmapFont font = GLHelvetica10BitmapFont.getInstance(gl);
            int numChars = font.getNumChars();
            for (int i = 0; i < nLabels; i++, md += labelInterval)
            {
                point = getPointAtMDepth((float) (i * labelInterval), true);
                float[] xyz = point.getXYZorT();
                if ((md % (5.0f * labelInterval)) != 0.0f)
                {
                    StsGLDraw.drawPoint(xyz, null, glPanel3d, 5, 1, 0.0f);
                }
                else
                {
                    StsGLDraw.drawPoint(xyz, null, glPanel3d, 10, 2, 0.0f);
                    float value = 0.0f;
                    if (drawLabelString.equals(MDEPTH_LABEL))
                    {
                        value = md;
                        //                       label = Float.toString(md);
                    }
                    else if (drawLabelString.equals(DEPTH_LABEL))
                    {
                        value = point.getZ();
                        //                        label = Float.toString(point.getZ());
                    }
                    else if (drawLabelString.equals(TIME_LABEL))
                    {
                        value = point.getT();
                        //                       label = Float.toString(point.getT());
                    }
                    label = labelFormat.format(value);
                    StsGLDraw.fontOutput(gl, xyz, label, font);
                }
            }
            glPanel3d.resetViewShift(gl);
        }
        displaySeismicCurtain(glPanel3d);
    }

    //TODO this is called on db reload, so we don't want to call dbFieldChanged every time; persist as a property?
    public void setLogTypeDisplay3dLeft(String logTypeName)
    {
         leftDisplayLogCurves = getAllLogCurvesOfType(logTypeName);
         // dbFieldChanged("leftDisplayLogCurves", leftDisplayLogCurves);
    }

    //TODO this is called on db reload, so we don't want to call dbFieldChanged every time; persist as a property?
    public void setLogTypeDisplay3dRight(String logTypeName)
    {
         rightDisplayLogCurves = getAllLogCurvesOfType(logTypeName);
    //     dbFieldChanged("rightDisplayLogCurves", rightDisplayLogCurves);
    }

    public void logDisplay3dLeft()
    {
        if(logCurves == null)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING,"There are no logs available for well " + getName());
            return;
        }
        StsObject[] logCurveObjects = logCurves.getElements();
        Object[] selectedLogCurveObjects = StsListSelectionDialog.getMultiSelectFromListDialog(currentModel.win3d,  "Log Selection", "Select left side log(s)", logCurveObjects);
        if(selectedLogCurveObjects.length == 0) return;
        leftDisplayLogCurves = (StsLogCurve[])StsMath.arraycopy(selectedLogCurveObjects, StsLogCurve.class);
        dbFieldChanged("leftDisplayLogCurves", leftDisplayLogCurves);
        currentModel.viewObjectRepaint(this, this);
    }

    public void logDisplay3dRight()
    {
        if(logCurves == null)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING,"There are no logs available for well " + getName());
            return;
        }
        StsObject[] logCurveObjects = logCurves.getElements();
        Object[] selectedLogCurveObjects = StsListSelectionDialog.getMultiSelectFromListDialog(currentModel.win3d,  "Log Selection", "Select right side log(s)", logCurveObjects);
        if(selectedLogCurveObjects == null) return;
        if(selectedLogCurveObjects.length == 0) return;
        rightDisplayLogCurves = (StsLogCurve[])StsMath.arraycopy(selectedLogCurveObjects, StsLogCurve.class);
        dbFieldChanged("rightDisplayLogCurves", rightDisplayLogCurves);
        currentModel.viewObjectRepaint(this, this);
    }

    public void logDisplay3d(StsLogCurve[] logCurves, float origin, StsGLPanel3d glPanel3d)
    {
        if (logCurves == null)return;
        for(StsLogCurve logCurve : logCurves)
            logCurve.display3d(glPanel3d, this, origin);
    }

    public void displaySeismicCurtain(StsGLPanel3d glPanel3d)
    {
        if(!getDrawCurtain()) return;
        if ((seismicCurtain != null) && (glPanel3d != null))
        {
            seismicCurtain.displayTexture3d(glPanel3d, drawCurtainTransparent);
            seismicCurtain.display(glPanel3d);
        }
    }

    public void createDisplaySeismicCurtain(StsGLPanel3d glPanel3d)
    {
        if(seismicCurtain == null)
        {
            StsSeismicVolume seismicVolume = (StsSeismicVolume)currentModel.getCurrentObject(StsSeismicVolume.class);
            if(seismicVolume == null) return;
            seismicCurtain = this.getCreateSeismicCurtain(seismicVolume);
        }
        displaySeismicCurtain(glPanel3d);
    }

    /** remove a well from the instance list and in the 3d window */
    public boolean delete()
    {
        // Remove well from any platforms it is assigned to
        StsPlatformClass pc = (StsPlatformClass) currentModel.getStsClass(StsPlatform.class);
        if (pc != null)
            pc.deleteWellFromPlatform(getName());

        if (!super.delete())
        {
            return false;
        }
        StsObjectRefList.deleteAll(logCurves);
        StsObjectRefList.deleteAll(lithZones);
        StsObjectRefList.deleteAll(markers);

        StsObject editTdSet = currentModel.getObjectWithName(StsEditTdSet.class, getName());
        if (editTdSet != null)
        {
            editTdSet.delete();
        }
        getWellClass().delete(this);
        return true;
    }

    public void deleteMarker(StsWellMarker wellMarker)
    {
        markers.delete(wellMarker);
    }

    public boolean addToProject()
    {
        StsProject project = currentModel.getProject();
        boolean hasTDCurve = getLastLogCurveOfType(StsWellKeywordIO.TIME) != null;
        if (zDomainSupported == TD_TIME_DEPTH && project.checkSetZDomain(TD_TIME_DEPTH, TD_DEPTH))
        {
            return super.addToProject(TD_TIME_DEPTH);
        }
        else if (project.checkSetZDomain(TD_DEPTH))
        {
            return super.addToProject(TD_DEPTH);
        }
        else
        {
            return false;
        }
    }

    // Log Curve routines

    public boolean constructWellDevCurves(StsLogVector[] logVectors, float curveNullValue, StsLogCurve tdCurve)
    {
        StsLogVector xVector = null, yVector = null, depthVector = null, mdepthVector = null, timeVector = null;
        float[] xArray, yArray, zArray, mArray, tArray;

        try
        {
            if (logVectors == null)
            {
                return false;
            }
            int nLogVectors = logVectors.length;

            StsLogVector[] logCurveVectors = new StsLogVector[0];
            for (int n = 0; n < nLogVectors; n++)
            {
                byte logVectorType = logVectors[n].type;
                if (logVectorType == StsLogVector.MDEPTH)
                {
                    mdepthVector = logVectors[n];
                    hasMDepths = true;
                }
                else if (logVectorType == StsLogVector.X)
                {
                    xVector = logVectors[n];
                }
                else if (logVectorType == StsLogVector.Y)
                {
                    yVector = logVectors[n];
                }
                else if (logVectorType == StsLogVector.DEPTH)
                {
                    depthVector = logVectors[n];
                }
                else if (logVectorType == StsLogVector.TIME)
                {
                    timeVector = logVectors[n];
                }
            }

            if (depthVector == null)
            {
                return false;
            }
            else
            {
                // Add datum Shift
                ;
            }

            if (mdepthVector == null)
            {
                mdepthVector = getMDepthsFromDepths(xVector, yVector, depthVector, StsLogVector.WELL_DEV_PREFIX,
                    depthVector.getUnits(), xVector.getUnits());
                if (mdepthVector == null)
                {
                    return false;
                }
            }
            else
            {
                // Add datum Shift
                ;
            }

            //		if(!setVectors(logVectors)) return false;
            //		StsLogCurve depthLogCurve = new StsLogCurve(well, mdepthVector, depthVector, null, 0);
            //		well.setDepthLogCurve(depthLogCurve);

            if (depthVector != null)
            {
                depthVector.checkMonotonic();
            }
            xArray = xVector.getFloats();
            if (xArray == null)
            {
                return false;
            }

            yArray = yVector.getFloats();
            if (yArray == null)
            {
                return false;
            }

            zArray = depthVector.getFloats();
            if (zArray == null)
            {
                return false;
            }

            mArray = null;
            if (mdepthVector != null)
            {
                mArray = mdepthVector.getFloats();
            }
            tArray = null;
            if (timeVector != null)
            {
                tArray = timeVector.getFloats();
            }

            int nValues = xArray.length;

            nativeHorizontalUnits = xVector.getUnits();
            nativeVerticalUnits = depthVector.getUnits();
            float hScalar = currentModel.getProject().getXyScalar(xVector.getUnits());
            float vScalar = currentModel.getProject().getDepthScalar(depthVector.getUnits());

            xOrigin = xVector.getOrigin() * hScalar;
            yOrigin = yVector.getOrigin() * hScalar;

            if (lineVertices == null)
            {
                lineVertices = StsObjectRefList.constructor(nValues, 1, "lineVertices", this);
            }

            if (timeVector != null)
            {
                StsPoint point = new StsPoint(5);
                for (int i = 0; i < nValues; i++)
                {
                    point.setX(xArray[i] * hScalar);
                    point.setY(yArray[i] * hScalar);
                    point.setZ(zArray[i] * vScalar);
                    point.setM(mArray[i] * vScalar);
                    point.setT(tArray[i]);
                    lineVertices.add(new StsSurfaceVertex(point, this));
                }
                setZDomainSupported(TD_TIME_DEPTH);
            }
            else if (currentModel.getProject().getSeismicVelocityModel() != null)
            {
                StsSeismicVelocityModel velocityModel = currentModel.getProject().getSeismicVelocityModel();
                //                velocityModel.checkInitializeIntervalVelocities();
                StsPoint point = new StsPoint(5);
                for (int i = 0; i < nValues; i++)
                {
                    point.setX(xArray[i] * hScalar);
                    point.setY(yArray[i] * hScalar);
                    point.setZ(zArray[i] * vScalar);
                    point.setM(mArray[i] * vScalar);
                    float t = (float) velocityModel.getT(point.v);
                    if (t != nullValue)
                    {
                        point.setT(t);
                    }
                    lineVertices.add(new StsSurfaceVertex(point, this));
                }
                setZDomainSupported(TD_TIME_DEPTH);
            }
            else if (tdCurve != null)
            {
                tArray = new float[nValues];
                for (int n = 0; n < nValues; n++)
                {
                    tArray[n] = tdCurve.getInterpolatedValue(zArray[n]);
                }

                StsPoint point = new StsPoint(5);
                for (int i = 0; i < nValues; i++)
                {
                    point.setX(xArray[i] * hScalar);
                    point.setY(yArray[i] * hScalar);
                    point.setZ(zArray[i] * vScalar);
                    point.setM(mArray[i] * vScalar);
                    point.setT(tArray[i]);
                    lineVertices.add(new StsSurfaceVertex(point, this));
                }
                setZDomainSupported(TD_TIME_DEPTH);
            }
            else
            {
                StsPoint point = new StsPoint(4);
                for (int i = 0; i < nValues; i++)
                {
                    point.setX(xArray[i] * hScalar);
                    point.setY(yArray[i] * hScalar);
                    point.setZ(zArray[i] * vScalar);
                    point.setM(mArray[i] * vScalar);
                    lineVertices.add(new StsSurfaceVertex(point, this));
                }
                setZDomainSupported(TD_DEPTH);
            }
            /** If we don't have a color yet, get it now */
            if (stsColor == null)
            {
                stsColor = currentModel.getCurrentSpectrumColor("Basic");
            }
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Unable to build well line geometry");
            return false;
        }
        numberOfElements = mdepthVector.getSize();
        return true;
    }

    public void setZDomainSupported(byte zDomain)
    {
        zDomainSupported = zDomain;
    }

    public void checkSetZDomainSupported(byte zDomain)
    {
        if (zDomainSupported == StsProject.TD_TIME_DEPTH) return;
        zDomainSupported = zDomain;
        dbFieldChanged("zDomainSupported", zDomain);
    }

    public boolean isInTime()
    {
        return zDomainSupported == TD_TIME_DEPTH;
    }
/*
    public void saveVertexTimesToDB(float[] times)
    {
        currentModel.addMethodCmd(this, "updateVertexTimes", new Object[]{times});
    }

    public void updateVertexTimes(float[] times)
    {
        int nVertices = lineVertices.getSize();
        for (int n = 0; n < nVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
            StsPoint point = vertex.getPoint();
            point.setT(times[n]);
        }
        computeXYZPoints();
    }
*/
    public void computeMarkerTimesFromMDepth(StsSeismicVelocityModel velocityModel)
    {
        if (markers == null) return;
        StsProject project = currentModel.getProject();
        int nMarkers = markers.getSize();
        for (int n = 0; n < nMarkers; n++)
        {
            StsWellMarker wellMarker = (StsWellMarker) markers.getElement(n);
            computeMarkerTimeFromMDepth(wellMarker, velocityModel);
        }
    }

    private boolean computeMarkerTimeFromMDepth(StsWellMarker wellMarker, StsSeismicVelocityModel velocityModel)
    {
        try
        {
            StsPoint location = wellMarker.getLocation();
            location.setT(getTimeFromMDepth(location.getM()));
            // wellMarker.dbFieldChanged("location", location);
            return true;
        }
        catch (Exception e)
        {
            return false;
        }
    }

    public boolean constructLogCurvesCheckVersions(StsLogVector[] logVectors, float curveNullValue)
    {
        try
        {
            if (logVectors == null)
            {
                return false;
            }
            int nLogVectors = logVectors.length;
            if (nLogVectors == 0)
            {
                return false;
            }

            // sort the logVectors by version
            ArrayList objects = new ArrayList(nLogVectors);
            for (int n = 0; n < nLogVectors; n++)
            {
                objects.add(logVectors[n]);

            }
            Comparator comparator = new VersionComparator();
            Collections.sort(objects, comparator);

            Iterator iter = objects.iterator();
            int currentVersion = StsParameters.nullInteger;
            int nVersionVectors = 0;
            while (iter.hasNext())
            {
                StsLogVector logVector = (StsLogVector) iter.next();
                int version = logVector.getVersion();
                if (version == currentVersion)
                {
                    logVectors[nVersionVectors++] = logVector;
                }
                else
                {
                    if (nVersionVectors > 0)
                    {
                        logVectors = (StsLogVector[]) StsMath.trimArray(logVectors, nVersionVectors);
                        StsLogCurve[] logCurves = StsLogCurve.constructLogCurves(this, logVectors, curveNullValue, version);
                        addLogCurves(logCurves);
                        nVersionVectors = 0;
                    }
                    currentVersion = version;
                    logVectors = new StsLogVector[nLogVectors];
                    logVectors[nVersionVectors++] = logVector;
                }
            }
            if (nVersionVectors > 0)
            {
                logVectors = (StsLogVector[]) StsMath.trimArray(logVectors, nVersionVectors);
                int version = logVectors[0].getVersion();
                StsLogCurve[] logCurves = StsLogCurve.constructLogCurves(this, logVectors, curveNullValue, version);
                addLogCurves(logCurves);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsLogCurve.constructWellDevCurvesCheckVersions() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    static public final class VersionComparator implements Comparator
    {
        VersionComparator()
        {
        }

        // order by versions and then order alphabetically
        public int compare(Object o1, Object o2)
        {
            StsLogVector v1 = (StsLogVector) o1;
            if (v1 == null)
            {
                return -1;
            }
            StsLogVector v2 = (StsLogVector) o2;
            if (v2 == null)
            {
                return 1;
            }

            int vv1 = v1.getVersion();
            int vv2 = v2.getVersion();

            // compa
            if (vv1 > vv2)
            {
                return 1;
            }
            if (vv1 < vv2)
            {
                return -1;
            }

            String s1 = v1.getName();
            String s2 = v2.getName();
            return s1.compareTo(s2);
        }
    }

    /** add to reference lists */

    public void addLogCurves(StsLogCurve[] logCurves)
    {
        for (int n = 0; n < logCurves.length; n++)
        {
            addLogCurve(logCurves[n]);
        }
    }

    public void addLogCurve(StsLogCurve logCurve)
    {
        if (logCurve == null) return;
        if (logCurves == null)
            logCurves = StsObjectRefList.constructor(10, 1, "logCurves", this);
        logCurve.setWell(this);
//         logCurve.checkAddMDepth(this);
        logCurves.add(logCurve);
    }

    public StsObjectRefList getLogCurves()
    {
        if (logCurves == null)
        {
            logCurves = StsObjectRefList.constructor(2, 2, "logCurves", this);
        }
        return logCurves;
    }

    public void deleteLogCurve(StsLogCurve logCurve)
    {
        logCurves.delete(logCurve);
        logCurve.setWell(null);
    }

    public StsLogCurve[] getUnusedLogCurves(StsLogCurve[] usedLogCurves)
    {
        StsLogCurve[] unusedLogCurves = null;

        int nLogCurves = logCurves.getSize();
        if (usedLogCurves == null || usedLogCurves.length == 0)
        {
            unusedLogCurves = new StsLogCurve[nLogCurves];
            for (int n = 0; n < nLogCurves; n++)
            {
                unusedLogCurves[n] = (StsLogCurve) logCurves.getElement(n);
            }
            return unusedLogCurves;
        }

        int nUsedLogCurves = usedLogCurves.length;
        int nUnusedLogCurves = nLogCurves - nUsedLogCurves;
        unusedLogCurves = new StsLogCurve[nUnusedLogCurves];
        int nn = 0;
        for (int n = 0; n < nLogCurves; n++)
        {
            StsLogCurve logCurve = (StsLogCurve) logCurves.getElement(n);
            if (logCurve.isInList(usedLogCurves))
            {
                continue;
            }
            unusedLogCurves[nn++] = logCurve;
        }
        return unusedLogCurves;
    }

    public int getNLogCurves()
    {
        return (logCurves == null) ? 0 : logCurves.getSize();
    }

    public StsLogCurve[] copyLogCurveArray()
    {
        if (logCurves == null)
            return new StsLogCurve[0];
        return (StsLogCurve[]) logCurves.copyArrayList(StsLogCurve.class);
    }

    /*
     int nCurves = getNLogCurves();
     if (nCurves==0) return new StsLogCurve[0];
     StsLogCurve[] curves = new StsLogCurve[nCurves];
     for (int i=0; i<nCurves; i++)
     {
      curves[i] = (StsLogCurve)logCurves.getElement(i);
     }
     return curves;
    }
    */
    public StsLogCurve[] getLogCurveList()
    {
        return copyLogCurveArray();
    }

    static public void printLogCurves(StsWell well) throws StsException
    {
        printLogCurves(well, 50);
    }

    /** print out log curve values */
    static public void printLogCurves(StsWell well, int increment)
    {
        StsLogCurve[] curves = null;
        try
        {
            curves = well.copyLogCurveArray();
        }
        catch (Exception e)
        {
            System.out.println("No log curves found.");
            return;
        }

        StringBuffer buffer = new StringBuffer("Well: " + well.getName() + "\n");
        buffer.append(" \n");
        if (curves == null)
        {
            buffer.append("No log curves found.\n");
        }
        else
        {
            String depthFmtPattern = "00000.00";
            String valueFmtPattern = "00000.0000";
            StsDecimalFormat depthFormat = new StsDecimalFormat(depthFmtPattern);
            StsDecimalFormat valueFormat = new StsDecimalFormat();
            for (int i = 0; i < curves.length; i++)
            {
                try
                {
                    String mdepthName = curves[i].getMDepthVector().getName();
                    String depthName = curves[i].getDepthVector().getName();
                    String curveName = curves[i].getName();
                    if (curveName.equals(mdepthName) || curveName.equals(depthName))
                        valueFormat.applyPattern(depthFmtPattern);
                    else
                        valueFormat.applyPattern(valueFmtPattern);
                    int nValues = curves[i].getValuesFloatVector().getSize();
                    buffer.append("Curve: " + curveName);
                    buffer.append("\nNumber of values: " + nValues + "\n \n");
                    buffer.append("index\t" + mdepthName + "\t\t" + depthName + "\t\t" + curveName + "\n");
                    StsFloatVector mdepthValues = curves[i].getMDepthFloatVector();
                    StsFloatVector depthValues = curves[i].getDepthFloatVector();
                    StsFloatVector curveValues = curves[i].getValuesFloatVector();
                    for (int j = 0; j < nValues; j += increment)
                    {
                        String mdepth = depthFormat.stripLeadingZeroes(mdepthValues.getElement(j));
                        String depth = depthFormat.stripLeadingZeroes(depthValues.getElement(j));
                        String value = valueFormat.stripLeadingZeroes(curveValues.getElement(j));
                        buffer.append(j + "\t" + mdepth + "\t" + depth + "\t" + value + "\n");
                    }
                    // print last value
                    if (nValues - 1 % increment != 0)
                    {
                        String mdepth = depthFormat.stripLeadingZeroes(mdepthValues.getElement(nValues - 1));
                        String depth = depthFormat.stripLeadingZeroes(depthValues.getElement(nValues - 1));
                        String value = valueFormat.stripLeadingZeroes(curveValues.getElement(nValues - 1));
                        buffer.append((nValues - 1) + "\t" + mdepth + "\t" + depth + "\t" + value + "\n");
                    }
                    buffer.append(" \n");
                }
                catch (Exception e)
                {
                    buffer.append("\nError in curve #" + i + ".  Continuing...\n");
                    continue;
                }
            }
        }

        // print out depth and curve names
        PrintWriter out = new PrintWriter(System.out, true); // needed for correct formatting
        out.println(buffer.toString());

        // display dialog box
        StsTextAreaDialog dialog = new StsTextAreaDialog(null, "Log Curve Listing for " + well.getName(),
            buffer.toString(), 40, 60, false);
        dialog.setLocationRelativeTo(getCurrentModel().win3d);
        dialog.setVisible(true);
    }

    public boolean hasLogCurves() { return logCurves != null && logCurves.getSize() > 0; }

    /** find a log curve in the list */
    public StsLogCurve getLastLogCurveOfType(String name)
    {
        if (name == null)
        {
            return null;
        }
        if (name.equalsIgnoreCase("none")) return null;
        if (logCurves == null)
        {
            return null;
        }
        //if(name.equalsIgnoreCase(StsWellClass.LOG_TYPE_NONE)) // Have to check here since "None" is the default alias
        //	return null;

        int nCurves = logCurves.getSize();
        for (int i = nCurves-1; i >= 0; i--)
        {
            StsLogCurve curve = (StsLogCurve) logCurves.getElement(i);
            if (curve.matchesName(name)) return curve;
        }
        return null;
    }

    public StsLogCurve[] getAllLogCurvesOfType(String curveTypeName)
    {
        if (curveTypeName == null)
        {
            return null;
        }
        if (curveTypeName.equalsIgnoreCase("none")) return null;
        if (logCurves == null)
        {
            return null;
        }
        //if(curveTypeName.equalsIgnoreCase(StsWellClass.LOG_TYPE_NONE)) // Have to check here since "None" is the default alias
        //	return null;

        int nCurves = logCurves.getSize();
        StsLogCurve[] matchingCurves = new StsLogCurve[0];
        for (int i = 0; i < nCurves; i++)
        {
            StsLogCurve curve = (StsLogCurve) logCurves.getElement(i);
            if (curve.logCurveTypeNameMatches(curveTypeName))
                matchingCurves = (StsLogCurve[])StsMath.arrayAddElement(matchingCurves, curve);

        }
        return matchingCurves;
    }

    public StsLogCurve getTdCurve()
    {
        return getLastLogCurveOfType(StsWellKeywordIO.TIME);
    }

    // Zone routines

    public void checkConstructWellZones()
    {
        if (zones != null)
        {
            return;
        }
        StsClass modelZones = currentModel.getCreateStsClass(StsZone.class); // these are zones in top to bottom order
        int nZones = zones.getSize();
        zones = new StsList(nZones, 1);
        for (int z = 0; z < nZones; z++)
        {
            StsZone zone = (StsZone) modelZones.getElement(z);
            StsWellZone wellZone = constructStratZone(zone);
            if (wellZone != null)
            {
                zones.add(wellZone);
            }
        }
    }

    private StsWellZone constructStratZone(StsZone zone)
    {
        StsModelSurface topHorizon = zone.getTopModelSurface();
        String topMarkerName = topHorizon.getMarkerName();
        StsWellMarker topMarker = getMarker(topMarkerName, STRAT);
        if (topMarker == null)
        {
            topMarker = StsWellMarker.constructor(topHorizon, this);
        }
        if (topMarker == null)
        {
            return null;
        }

        StsModelSurface baseHorizon = zone.getBaseModelSurface();
        String baseMarkerName = baseHorizon.getMarkerName();
        StsWellMarker baseMarker = getMarker(baseMarkerName, STRAT);
        if (baseMarker == null)
        {
            baseMarker = StsWellMarker.constructor(baseHorizon, this);
        }
        if (baseMarker == null)
        {
            return null;
        }

        return StsWellZone.constructor(this, STRAT, topMarkerName, topMarker, baseMarker, zone);
    }

    private StsWellMarker getMarker(String markerName, int type)
    {
        if (markers == null)
        {
            return null;
        }
        int nMarkers = markers.getSize();
        for (int n = 0; n < nMarkers; n++)
        {
            StsWellMarker marker = (StsWellMarker) markers.getElement(n);
            if (marker.getType() == type && marker.getName() == markerName)
            {
                return marker;
            }
        }
        return null;
    }

    public void addZone(StsWellZone zone)
    {
        if (zone.getZoneType() == StsWellZoneSet.STRAT)
        {
            if (zones == null)
            {
                zones = new StsList(10, 1);

            }
            zones.add(zone);
        }
        else if (zone.getZoneType() == StsWellZoneSet.LITH)
        {
            if (lithZones == null)
            {
                lithZones = StsObjectRefList.constructor(10, 1, "lithZones", this);

            }
            lithZones.add(zone);
        }
        //        setWellLineNeedsRebuild();
    }

    /* get an number of well zones (doesn't handle fault zones) */
    public int getNZones(int type)
    {
        if (type == StsWellZoneSet.STRAT && zones != null)
        {
            return zones.getSize();
        }
        else if (type == StsWellZoneSet.LITH && lithZones != null)
        {
            return lithZones.getSize();
        }
        else
        {
            return 0;
        }
    }

    /* get an array of well zones (doesn't handle fault zones) */
    public Object[] getZoneArray(int type)
    {
        if (type == StsWellZoneSet.STRAT && zones != null)
        {
            return zones.getList();
        }
        else if (type == StsWellZoneSet.LITH && lithZones != null)
        {
            return lithZones.getElements();
        }
        else
        {
            return null;
        }
    }

    /* get an list of well zones (doesn't handle fault zones) */
    public String[] getZoneList(int type)
    {
        Object[] zoneArray = getZoneArray(type);
        if (zoneArray == null)
        {
            return null;
        }
        String[] zoneList = new String[zoneArray.length];
        for (int i = 0; i < zoneArray.length; i++)
        {
            StsLineZone zone = (StsLineZone)zoneArray[i];
            zoneList[i] = zone.getLabel();
        }
        return zoneList;
    }

    public StsWellZone getZone(String name, int type)
    {
        if (name == null)
        {
            return null;
        }
        Object[] zoneArray = getZoneArray(type);
        if (zoneArray == null)
        {
            return null;
        }
        for (int i = 0; i < zoneArray.length; i++)
        {
            StsLineZone zone = (StsLineZone)zoneArray[i];
            if (name.equals(zone.getLabel()))
                return (StsWellZone) zoneArray[i];
        }
        return null;
    }
    
    public void setStsColor(StsColor color)
    {
        super.setStsColor(color);
        currentModel.viewObjectRepaint(this, this);
    }
    public float getDepthFromMDepth(float mdepth)
    {
        return StsMath.interpolateValue(rotatedPoints, mdepth, 3, 2);
    }

    public float getTimeFromMDepth(float mdepth)
    {
        StsSeismicVelocityModel velocityModel = currentModel.getProject().getSeismicVelocityModel();
        if(velocityModel != null)
        {
            
        }
        float depth = getDepthFromMDepth(mdepth);

        return StsMath.interpolateValue(rotatedPoints, mdepth, 3, 4);
    }

    public float getDepthFromTime(float mdepth)
    {
        return StsMath.interpolateValue(rotatedPoints, mdepth, 4, 2);
    }

    // Not reliable for horizontal well
    public StsPoint getPointFromDepth(float depth)
    {
        return StsMath.interpolatePoint(rotatedPoints, depth, 2);
    }

    public StsPoint getPointFromDepth(float depth, boolean extrapolate)
    {
        return StsMath.interpolatePoint(depth, rotatedPoints, 2, extrapolate);
    }

    public void setIsVisible(boolean b)
    {
        if (b == isVisible)
        {
            return;
        }
        isVisible = b;
        dbFieldChanged("isVisible", isVisible);
        currentModel.win3dDisplayAll();
    }

    /** get measured depth from tvd by linear interpolation of tvd */
    public float getMDepthFromDepth(float depth)
    {
        if (rotatedPoints != null)
        {
            return StsMath.interpolateValue(rotatedPoints, depth, 2, 3);
        }
        if (lineVertices == null || lineVertices.getSize() == 0) return depth;

        computePoints();
        if (rotatedPoints == null) return depth;
        return StsMath.interpolateValue(rotatedPoints, depth, 2, 3);
    }

    /** get measured depth from tvd by linear interpolation of tvd */
    public float getMDepthFromTime(float time)
    {
        if (rotatedPoints != null)
        {
            return StsMath.interpolateValue(rotatedPoints, time, 4, 3);
        }
        if (lineVertices == null)
        {
            return time;
        }
        int nLineVertices = lineVertices.getSize();
        if (nLineVertices == 0)
        {
            return time;
        }
        computePoints();
        if (rotatedPoints == null) return time;
        return StsMath.interpolateValue(rotatedPoints, time, 4, 3);
    }

    /** get interpolated value for one array given another */
    static public float getInterpolatedValue(float av, float[] a, float[] b)
    {
        if (a == null || b == null || a.length != b.length)
        {
            return nullValue;
        }
        return StsMath.interpolateValue(av, a, b, a.length);
    }

    /** convert a vector of measured depths to z depths */
    public StsLogVector getDepthsFromMDepths(StsLogVector mdepths, String binDir)
    {
        try
        {
            float[] mdepthValues = mdepths.getFloats();
            int nValues = mdepthValues.length;
            float[] depthValues = new float[nValues];
            for (int i = 0; i < nValues; i++)
            {
                depthValues[i] = getDepthFromMDepth(mdepthValues[i]);
            }
            String asciiFilename = null;
            String binaryFilename = new String("well-logs.bin." + getName() + ".DEPTH.0");
            StsLogVector depthVector = new StsLogVector(asciiFilename, binaryFilename, "DEPTH");
            depthVector.setValues(new StsFloatVector(depthValues));
            if (!depthVector.hasBinaryFile(binDir))
            {
                depthVector.checkWriteBinaryFile(binDir);
            }
            return depthVector;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWelLine.getDepthsFromMDepths() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsLogVector getDepthsFromMDepths(StsLogVector mdepths)
    {
        try
        {
            float[] mdepthValues = mdepths.getFloats();
            int nValues = mdepthValues.length;
            float[] depthValues = new float[nValues];
            for (int i = 0; i < nValues; i++)
                depthValues[i] = getDepthFromMDepth(mdepthValues[i]);
            StsLogVector depthVector = new StsLogVector(depthValues);
            return depthVector;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getDepthsFromMDepths", e);
            return null;
        }
    }

    public double[] getMDepthsFromTimes(float tMin, float tInc, int nValues)
    {
        try
        {
            double[] mdepthValues = new double[nValues];
            float time = tMin;
            for (int n = 0; n < nValues; n++, time += tInc)
                mdepthValues[n] = getMDepthFromTime(time);
            return mdepthValues;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWelLine.getMDepthsFromTimes() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public double[] getMDepthsFromDepths(float zMin, float zInc, int nValues)
    {
        try
        {
            double[] mdepthValues = new double[nValues];
            float z = zMin;
            for (int n = 0; n < nValues; n++, z += zInc)
                mdepthValues[n] = getMDepthFromDepth(z);
            return mdepthValues;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWelLine.getMDepthsFromTimes() failed.", e, StsException.WARNING);
            return null;
        }
    }

    /** convert a vector of z depths to measured depths */
    public StsLogVector getMDepthsFromDepths(StsLogVector xVector, StsLogVector yVector, StsLogVector zVector,
                                             String group, byte vUnits, byte hUnits)
    {
        try
        {
            if (xVector == null || yVector == null || zVector == null)
            {
                return null;
            }

            float[] xArray = xVector.getFloats();
            float[] yArray = yVector.getFloats();
            float[] zArray = zVector.getFloats();

            int nValues = xArray.length;
            float[] mdepthArray = new float[nValues];
            mdepthArray[0] = zArray[0];
            float x2 = xArray[0];
            float y2 = yArray[0];
            float z2 = zArray[0];

            for (int i = 1; i < nValues; i++)
            {
                float x1 = x2;
                x2 = xArray[i];
                float dx = x2 - x1;
                float y1 = y2;
                y2 = yArray[i];
                float dy = y2 - y1;
                float z1 = z2;
                z2 = zArray[i];
                float dz = z2 - z1;
                float dm = (float) Math.sqrt(dx * dx + dy * dy + dz * dz);
                mdepthArray[i] = mdepthArray[i - 1] + dm;
            }

            String asciiFilename = null;
            String binaryFilename = new String(group + ".bin." + getName() + ".MDEPTH.0");
            StsLogVector mdepthVector = new StsLogVector(asciiFilename, binaryFilename, "MDEPTH");
            mdepthVector.setValues(new StsFloatVector(mdepthArray));
            mdepthVector.setUnits(vUnits);
            mdepthVector.checkWriteBinaryFile(currentModel.getProject().getBinaryFullDirString());
            return mdepthVector;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellLine.getMDepthsFromDepths() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsPoint getPointFromLocation(float x, float y, float depth)
    {
        StsPoint point = StsMath.getNearestPointOnLine(new StsPoint(x, y, depth), rotatedPoints, 3);
        return point;
    }

    public StsLogVector getDepthsFromMDepthsUsingWellDev(float[] mdepths)
    {
        if (mdepths == null || lineVertices == null)
        {
            return null;
        }
        float[] depthArray = getDepthArrayFromVertices();
        float[] mdepthArray = getMDepthArrayFromVertices();
        int nDepths = mdepths.length;
        float[] depths = new float[nDepths];
        for (int n = 0; n < nDepths; n++)
        {
            depths[n] = StsMath.interpolateValue(mdepths[n], mdepthArray, depthArray);
        }
        return new StsLogVector(StsLogVector.DEPTH, depths);
    }

    public StsLogVector getMDepthsVectorFromDepths(float[] depths)
    {
        if (depths == null || lineVertices == null)
        {
            return null;
        }
        float[] depthArray = getDepthArrayFromVertices();
        float[] mdepthArray = getMDepthArrayFromVertices();
        int nDepths = depths.length;
        float[] mdepths = new float[nDepths];
        for (int n = 0; n < nDepths; n++)
        {
            mdepths[n] = StsMath.interpolateValue(depths[n], depthArray, mdepthArray);
        }
        return new StsLogVector(StsLogVector.MDEPTH, mdepths);
    }

    public float[] getMDepthsFromDepths(float[] depths)
    {
        if (depths == null || lineVertices == null) return null;
        float[] depthArray = getDepthArrayFromVertices();
        float[] mdepthArray = getMDepthArrayFromVertices();
        int nDepths = depths.length;
        float[] mdepths = new float[nDepths];
        for (int n = 0; n < nDepths; n++)
        {
            mdepths[n] = StsMath.interpolateValue(depths[n], depthArray, mdepthArray);
        }
        return mdepths;
    }

    public boolean checkAddMDepthToDev(StsLogVector[] curveLogVectors)
    {
        StsLogVector mDepthVector = StsLogVector.getVectorOfType(curveLogVectors, StsLogVector.MDEPTH);
        StsLogVector depthVector = StsLogVector.getVectorOfType(curveLogVectors, StsLogVector.DEPTH);
        if (mDepthVector == null || depthVector == null)
        {
            return false;
        }
        return checkAddMDepthVector(mDepthVector, depthVector);
    }

    /**
     * Measured depths aren't available from dev survey, so use this log curves which
     * generally have measured depths.  However, logged interval doesn't cover the entire
     * well, so adjust measured depths above and below by offseting the estimated measured
     * depths computed when the dev survey was read so that they tie at the logged interval.
     */

    public boolean checkAddMDepthVector(StsLogVector mDepthVector, StsLogVector depthVector)
    {
        int nLineVertices = lineVertices.getSize();
        StsPoint[] points = new StsPoint[nLineVertices];
        for (int n = 0; n < nLineVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
            points[n] = vertex.getPoint();

        }
        float vScalar = currentModel.getProject().getDepthScalar(mDepthVector.getUnits());
        float[] curveMDepths = mDepthVector.getFloats();
        float[] curveDepths = depthVector.getFloats();
        for (int j = 0; j < curveMDepths.length; j++)
        {
            curveMDepths[j] = curveMDepths[j] * vScalar;
            curveDepths[j] = curveDepths[j] * vScalar;
        }
        float[] newMDepths = new float[nLineVertices];
        float mDepthError = 0.0f;

        // get the measured offset adjustment above and below the logged interval
        float curveDepthMin = curveDepths[0];
        float devMDepth = getMDepthFromDepth(curveDepthMin);
        float aboveOffset = (curveMDepths[0]) - devMDepth;
        int n = 0;
        // adjust dev mDepths above the first logCurve depth
        for (; n < nLineVertices; n++)
        {
            float devDepth = points[n].getZ();
            if (devDepth >= curveDepthMin)
            {
                break;
            }
            devMDepth = points[n].getM();
            newMDepths[n] = devMDepth + aboveOffset;
        }
        int nLastCurveValue = curveMDepths.length - 1;
        float curveDepthMax = curveDepths[nLastCurveValue];
        for (; n < nLineVertices; n++)
        {
            float devDepth = points[n].getZ();
            if (devDepth > curveDepthMax)
            {
                break;
            }
            newMDepths[n] = StsMath.interpolateValue(points[n].getZ(), curveDepths, curveMDepths);
            if (hasMDepths)
            {
                mDepthError = Math.max(mDepthError, (Math.abs(newMDepths[n] - points[n].getM())));
            }
        }
        devMDepth = getMDepthFromDepth(curveDepthMax);
        float belowOffset = curveMDepths[nLastCurveValue] - devMDepth;
        for (; n < nLineVertices; n++)
        {
            devMDepth = points[n].getM();
            newMDepths[n] = devMDepth + belowOffset;
        }
        if (mDepthError > maxMDepthError)
        {
            boolean answer = StsYesNoDialog.questionValue(currentModel.win3d, "Measured depth values differ between " +
                "deviation survey files and log curve files: error " + mDepthError + " max allowed " + maxMDepthError);
            if (answer == false)
            {
                return false;
            }
        }
        for (n = 0; n < nLineVertices; n++)
        {
            points[n].setM(newMDepths[n]);
        }

        hasMDepths = true;
        return true;
    }

    /* get an array of well zones that need display (doesn't handle fault zones) */
    public Iterator getDisplayZoneArray(int type)
    {
        if (type == StsWellZoneSet.STRAT && zones != null)
        {
            return new DisplayZoneIterator(zones);
        }
        else if (type == StsWellZoneSet.LITH && lithZones != null)
        {
            return new DisplayZoneIterator(lithZones.getList());
        }
        return new DisplayZoneIterator();
    }

    class DisplayZoneIterator implements Iterator
    {
        StsList zones;
        StsWellZone next = null;
        int n = 0;

        DisplayZoneIterator()
        {
        }

        DisplayZoneIterator(StsList zones)
        {
            this.zones = zones;
            setNext();
        }

        private void setNext()
        {
            if (zones == null)
            {
                return;
            }
        }

        public boolean hasNext()
        {
            return next != null;
        }

        public Object next()
        {
            StsWellZone current = next;
            while ((next = (StsWellZone) zones.getElement(n++)) != null)
            {
                if (next.getDisplayOnWellLine())
                {
                    break;
                }
            }
            return current;
        }

        public void remove()
        {}
    }

    public void setDrawMarkers(boolean b)
    {
        if (drawMarkers == b)
        {
            return;
        }
        drawMarkers = b;
        currentModel.win3dDisplay();
    }

    public boolean getDrawMarkers()
    {
        return drawMarkers;
    }

    public void setDrawPerfMarkers(boolean b)
    {
        if (drawPerfMarkers == b)
        {
            return;
        }
        drawPerfMarkers = b;
        currentModel.win3dDisplay();
    }

    public boolean getDrawPerfMarkers()
    {
        return drawPerfMarkers;
    }

    public void setDrawEquipmentMarkers(boolean b)
    {
        if (drawEquipmentMarkers == b)
        {
            return;
        }
        drawEquipmentMarkers = b;
        currentModel.win3dDisplay();
    }

    public boolean getDrawEquipmentMarkers()
    {
        return drawEquipmentMarkers;
    }

    public boolean getDrawCurtainTransparent()
    {
        return drawCurtainTransparent;
    }

    public void setDrawCurtainTransparent(boolean value)
    {
        drawCurtainTransparent = value;
        dbFieldChanged("drawCurtainTransparent", drawCurtainTransparent);
    }

    public boolean getDrawCurtain()
    {
        return isDrawingCurtain;
    }

    public void setDrawCurtain(boolean curtain)
    {
        if (isDrawingCurtain == curtain)
        {
            return;
        }
        isDrawingCurtain = curtain;
        dbFieldChanged("isDrawingCurtain", isDrawingCurtain);
        if (isDrawingCurtain)
            createCurtain();
        else
        {
            deleteSeismicCurtain();
            currentModel.win3dDisplay();
        }
    }

    public void createCurtain()
    {
        Runnable runCreateCurtain = new Runnable()
        {
            public void run()
            {
                createSeismicCurtain();
                currentModel.win3dDisplay();
            }
        };
        Thread runCreateCurtainThread = new Thread(runCreateCurtain);
        runCreateCurtainThread.start();
    }

    public String getDrawLabelString()
    {
        return drawLabelString;
    }

    public void setDrawLabelString(String labelString)
    {
        drawLabelString = labelString;
        dbFieldChanged("drawLabelString", drawLabelString);
        currentModel.win3dDisplay();
        return;
    }

    public float getLabelInterval()
    {
        return labelInterval;
    }

    public void setLabelInterval(float value)
    {
        labelInterval = value;
        dbFieldChanged("labelInterval", labelInterval);
        currentModel.win3dDisplay();
        return;
    }

    /** marker methods */

    public void addModelSurfaceMarkers()
    {
        StsObject[] surfaces = currentModel.getObjectList(StsModelSurface.class);
        int nSurfaces = surfaces.length;
        for (int n = 0; n < nSurfaces; n++)
        {
            StsModelSurface surface = (StsModelSurface) surfaces[n];
            if (getMarker(surface.getName()) != null)
            {
                continue; // already have it
            }
            try
            {
                StsWellMarker.constructor(surface, this);
            }
            catch (Exception e)
            {}
        }
    }

    public void addMarker(StsWellMarker marker)
    {
        if (markers == null)
        {
            markers = StsObjectRefList.constructor(10, 1, "markers", this);
        }
        markers.add(marker);
    }

    public StsObjectRefList getMarkers()
    {
        return markers;
    }

    public boolean hasMarkers()
    {
        return getNMarkers() > 0;
    }

    public int getNMarkers()
    {
        return (markers == null) ? 0 : markers.getSize();
    }

    public StsWellMarker[] getMarkerArray()
    {
        int nMarkers = getNMarkers();
        if (nMarkers == 0)
        {
            return null;
        }
        StsWellMarker[] markerArray = new StsWellMarker[nMarkers];
        for (int i = 0; i < nMarkers; i++)
        {
            markerArray[i] = (StsWellMarker) markers.getElement(i);
        }
        return markerArray;
    }

    public String[] getMarkerList()
    {
        StsWellMarker[] markerArray = getMarkerArray();
        if (markerArray == null)
        {
            return null;
        }
        String[] markerList = new String[markerArray.length];
        for (int i = 0; i < markerArray.length; i++)
        {
            markerList[i] = markerArray[i].getName();
        }
        return markerList;
    }

    public StsWellMarker getMarker(String name)
    {
        if (name == null)
        {
            return null;
        }
        int nMarkers = getNMarkers();
        if (nMarkers == 0)
        {
            return null;
        }
        for (int i = 0; i < nMarkers; i++)
        {
            StsWellMarker marker = (StsWellMarker) markers.getElement(i);
            if (name.equals(marker.getName()))
            {
                return marker;
            }
        }
        return null;
    }

    public void adjustTimes(float[] adjustedTimes)
    {
        if (adjustedTimes == null)
        {
            return;
        }
        StsPoint[] lineVertexPoints = getLineVertexPoints();
        for (int n = 0; n < lineVertexPoints.length; n++)
        {
            lineVertexPoints[n].setT(adjustedTimes[n]);
        }
        computePoints();
        adjustMarkerTimes();
    }

    public void adjustMarkerTimes()
    {
        if (markers == null)
        {
            return;
        }
        int nMarkers = markers.getSize();
        for (int n = 0; n < nMarkers; n++)
        {
            StsWellMarker marker = (StsWellMarker) markers.getElement(n);
            marker.adjustTime();
        }
    }

    /** Build a wellViewModel and a frameView if we don't have one. If we have a wellViewModel, but no frameView, built that. */
    public void openOrPopWindow()
    {
        if(wellViewModel == null)
            wellViewModel = new StsWellFrameViewModel(this);
        else if(wellViewModel.wellWindowFrame == null)
            wellViewModel.buildFrameView();
    }

    public void close()
    {
        wellViewModel.wellWindowFrame = null;
    }

    public void setWellFrameViewModel(StsWellFrameViewModel wellViewModel)
    {
        this.wellViewModel = wellViewModel;
    }

    public void pick(GL gl, StsGLPanel glPanel)
    {
        super.pick(gl, glPanel);
    }

    public void mouseSelectedEdit(StsMouse mouse)
    {
        logMessage();
        openOrPopWindow();
    }

    public boolean initialize(StsModel model)
    {
        return initialize();
    }

    /**
     * Initialize well even if on uninitialized section. Return true only if
     * not on section or section is initialized
     */
    public boolean initialize()
    {
        if (initialized)return true;
        if (!super.initialize())return false;
        if (isDrawingCurtain) getCreateSeismicCurtain();
        return true;
    }

    public void logMessage()
    {
        logMessage("Well: " + this.getName() + " " + lineOnSectionLabel());
    }

    /*
    public float[] getVerticesTimes()
    {
     if(verticesTimes != null) return verticesTimes;
     StsLogCurve timeCurve = getLogCurve(TIME);
     float[] wellDepths = getVerticesDepths();
     int nValues = wellDepths.length;
     verticesTimes = new float[nValues];
     for(int n = 0; n < nValues; n++)
      verticesTimes[n] = timeCurve.interpolatedValue(wellDepths[n]);
     return verticesTimes;
    }

    public float[] getVerticesDepths()
    {
     if(verticesDepths != null) return verticesDepths;
     verticesDepths = this.getDepthArrayFromVertices();
     return verticesDepths;
    }

    public float[] getVerticesMDepths()
    {
     if(verticesMDepths != null) return verticesMDepths;
     verticesMDepths = getMDepthArrayFromVertices();
     return verticesDepths;
    }
    */
    /*
    public StsPoint[] getTimePoints()
    {
     if(timePoints != null) return timePoints;

     StsLogCurve timeCurve = getLogCurve(TIME);
     if(timeCurve == null) return null;
     float[] depths = getVerticesDepths();
     float[] times = getVerticesTimes();
     int nPoints = points.length;
     timePoints = new StsPoint[nPoints];
     for(int n = 0; n < nPoints; n++)
     {
      timePoints[n] = new StsPoint(depthPoints[n]);
      float time = StsLogCurve.interpolatedValue(depths, times, depthPoints[n].getZ());
      timePoints[n].setZ(time);
     }
     return timePoints;
    }
    */
    public float getValueFromMDepth(float mdepth, String typeString)
    {
        if (typeString == StsWellViewModel.MDEPTH)
        {
            return mdepth;
        }
        if (typeString == StsWellViewModel.DEPTH)
        {
            return getDepthFromMDepth(mdepth);
        }
        else if (typeString == StsWellViewModel.TIME)
        {
            return getTimeFromMDepth(mdepth);
        }
        else
        {
            return nullValue;
        }
    }

    public float getMinMDepth()
    {
        StsPoint point = getPoint(0);
        if (point == null)
        {
            return nullValue;
        }
        return point.getM();
    }

    public StsPoint getPoint(int index)
    {
        if (rotatedPoints != null && index < rotatedPoints.length)
        {
            return rotatedPoints[index];
        }
        else if (lineVertices != null && index < lineVertices.getSize())
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(index);
            return vertex.getPoint();
        }
        return null;
    }

    public StsPoint getLastPoint()
    {
        StsPoint[] points = getRotatedPoints();
        if (points != null)
            return points[points.length - 1];
        else if (lineVertices != null && lineVertices.getSize() > 0)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getLast();
            return vertex.getPoint();
        }
        return null;
    }

    public float getMaxMDepth()
    {
        StsPoint point = getLastPoint();
        if (point == null)
        {
            return nullValue;
        }
        return point.getM();
    }

    public StsPoint[] getPointsFromLineVertices()
    {
        if (lineVertices == null)
        {
            return null;
        }
        int nLineVertices = lineVertices.getSize();
        if (nLineVertices == 0)
        {
            return null;
        }
        rotatedPoints = new StsPoint[nLineVertices];
        for (int n = 0; n < nLineVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
            rotatedPoints[n] = vertex.getPoint();
        }
        return rotatedPoints;
    }

    public void createSeismicCurtain()
    {
		if(currentModel == null || currentModel.win3d == null) return;
        StsObject object = ((StsComboBoxToolbar)currentModel.win3d.getToolbarNamed(StsComboBoxToolbar.NAME)).getSelectedObject();
        if(object instanceof StsVirtualVolume)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Seismic curtains only support seismic volume, not virtual volumes.\n\nVirtual volume can be exported and re-imported as a seismic volume.");
            setDrawCurtain(false);
            getObjectPanel().updateBeans();
        }
        else
            createSeismicCurtain(((StsSeismicVolume) object));
    }

    public void createSeismicCurtain(StsSeismicVolume vol)
    {
        if (seismicCurtain != null)
        {
            if (seismicCurtain.getSeismicVolume() == vol)
                return;
        }
        if(vol == null)
        {
            new StsMessage(currentModel.win3d, StsMessage.ERROR, "Well curtain requires a seismic volume to be loaded into the project.");
            isDrawingCurtain = false;
            return;
        }
        else
        {
            StsPoint[] rotatedPoints = getRotatedPoints();
            isDrawingCurtain = true;
            if (isVertical)
                seismicCurtain = new StsSeismicCurtain(currentModel, rotatedPoints, vol);
            else
            {
                // Pass the line for the curtain extraction
                StsPoint[] shiftedPoints = shiftPoints(rotatedPoints, curtainOffset);
                seismicCurtain = new StsSeismicCurtain(currentModel, shiftedPoints, vol);
            }
        }
    }

    /** shift points along normal to line between first and last points */
    private StsPoint[] shiftPoints(StsPoint[] rotatedPoints, float offset)
    {
        float[] normal = StsMath.horizontalNormal(rotatedPoints[0].getPointXYZ(), rotatedPoints[rotatedPoints.length - 1].getXYZ(), 1);
        float xShift = normal[0] * offset;
        float yShift = normal[1] * offset;
        StsPoint[] shiftedPoints = StsPoint.copy(rotatedPoints);
        for (int n = 0; n < shiftedPoints.length; n++)
        {
            shiftedPoints[n].v[0] += xShift;
            shiftedPoints[n].v[1] += yShift;
        }
        return shiftedPoints;
    }

    public void deleteSeismicCurtain()
    {
        if (seismicCurtain != null)
        {
            seismicCurtain.delete();
            //           currentModel.win3dDisplay();
        }
        seismicCurtain = null;
    }

	public void addTimeDepthToVels(StsSeismicVelocityModel velocityModel)
	{

    }

    public void checkAdjustWellTimes(StsSeismicVelocityModel velocityModel)
    {
        // domain should always be depth actually
        //        if (zDomainOriginal == StsParameters.TD_TIME)
        {
            adjustFromVelocityModel(velocityModel);
            adjustNonSurfaceWellMarkerTimes(velocityModel);
        }
        /*
        else
        {
            adjustWellTimes(velocityModel);
            convertMarkersToTime(velocityModel);
        }
        */
        checkSetZDomainSupported(StsParameters.TD_TIME_DEPTH);
    }

    /*
          public boolean checkConvertToTime(StsSeismicVelocityModel velocityModel)
         {
             if(lineVertices == null) return false;
             int nVertices = lineVertices.getSize();
             StsProject project = currentModel.getProject();
             float[] times = new float[nVertices];
             for (int n = 0; n < nVertices; n++)
             {
                 StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
                 StsPoint point = vertex.getPoint();
                StsPoint newPoint = adjustPointTime(point, project, velocityModel);
                if (newPoint == null)
                {
                    StsMessageFiles.errorMessage("Failed to convert well " + getName() + " to time.");
                    return false;
                }
                vertex.setPoint(newPoint);
                times[n] = newPoint.getT();
             }
             saveVertexTimesToDB(times);
     //       currentModel.addMethodCmd(this, "computeXYZPoints", new Object[0] );
             computeXYZPoints(); // generate splined points between vertices
             setZDomainSupported(TD_TIME_DEPTH);
             dbFieldChanged("zDomainSupported", TD_TIME_DEPTH);
     //        convertMarkersToTime(velocityModel);
             return true;
         }
    */
    public void saveVertexDepthsToDB(float[] depths)
    {
        currentModel.addMethodCmd(this, "updateVertexDepths", new Object[]
            {depths});
    }

    public void updateVertexDepths(float[] depths)
    {
        int nVertices = lineVertices.getSize();
        for (int n = 0; n < nVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
            StsPoint point = vertex.getPoint();
            point.setZ(depths[n]);
        }
        computePoints();
    }

    public void adjustWellDepths(StsSeismicVelocityModel velocityModel)
    {
        // get vertex points in rotated project coordinate system
        StsPoint[] linePoints = computeRotatedCoorVertexPoints();
        int nLineVertices = lineVertices.getSize();
        float maxChange = 0.0f;
        float[] depths = new float[nLineVertices];
        for (int v = 0; v < nLineVertices; v++)
        {
            StsSurfaceVertex lineVertex = (StsSurfaceVertex) lineVertices.getElement(v);
            StsPoint linePoint = lineVertex.getPoint();
            // compute t for rotated coordinate point and put it in the unrotated point
            float z;
            try
            {
                float[] xyzmt = linePoints[v].v;
                z = (float) velocityModel.getZ(xyzmt[0], xyzmt[1], xyzmt[4]);
            }
            catch (Exception e)
            {
                StsMessageFiles.errorMessage("Failed to adjust well " + getName() + " points probably not in time.");
                return;
            }
            float currentZ = linePoints[v].getZ();
            maxChange = Math.max(maxChange, Math.abs(z - currentZ));
            linePoint.setZ(z);
            depths[v] = z;
            //            StsGridSectionPoint surfacePoint = lineVertex.getSurfacePoint();
            //            surfacePoint.dbFieldChanged("point", linePoint);
        }
        saveVertexDepthsToDB(depths);
        //       currentModel.addMethodCmd(this, "computeXYZPoints", new Object[0]);
        StsMessageFiles.logMessage("Well " + getName() + " max depth adjustment " + maxChange);
        adjustDepthPoints(velocityModel);
        currentModel.addMethodCmd(this, "adjustWellDepthPoints", new Object[]
            {velocityModel}, "adjustWellDepthPoints for " + getName());
    }

    /*
        public void adjustWellDepths(StsSeismicVelocityModel velocityModel)
         {
             // get vertex points in rotated project coordinate system
             StsPoint[] linePoints = computeProjectPoints();
             int nLineVertices = lineVertices.getSize();
             float maxChange = 0.0f;
             float[] depths = new float[nLineVertices];
             for (int v = 0; v < nLineVertices; v++)
             {
                 StsSurfaceVertex lineVertex = (StsSurfaceVertex) lineVertices.getElement(v);
                 StsPoint linePoint = lineVertex.getPoint();
                 // compute t for rotated coordinate point and put it in the unrotated point
                 float z;
                 try
                 {
                     float[] xyzmt = linePoints[v].v;
                     z = (float) velocityModel.getZ(xyzmt[0], xyzmt[1], xyzmt[4]);
                 }
                 catch (Exception e)
                 {
                     StsMessageFiles.errorMessage("Failed to adjust well " + getName() + " points probably not in time.");
                     return;
                 }
                 float currentZ = linePoints[v].getZ();
                 maxChange = Math.max(maxChange, Math.abs(z - currentZ));
                 linePoint.setZ(z);
                 depths[v] = z;
    //            StsGridSectionPoint surfacePoint = lineVertex.getSurfacePoint();
    //            surfacePoint.dbFieldChanged("point", linePoint);
             }
             saveVertexDepthsToDB(depths);
             //       currentModel.addMethodCmd(this, "computeXYZPoints", new Object[0]);
             StsMessageFiles.logMessage("Well " + getName() + " max depth adjustment " + maxChange);
             adjustWellDepthPoints(velocityModel);
             currentModel.addMethodCmd(this, "adjustWellDepthPoints", new Object[]
                                       {velocityModel}, "adjustWellDepthPoints for " + getName());
        }
    */
    public void adjustNonSurfaceWellMarkerTimes(StsSeismicVelocityModel velocityModel)
    {
        float oldTime = 0.0f;
        float time;
        StsWellMarker wellMarker;

        StsObjectRefList markers = getMarkers();
        if (markers == null)
        {
            return;
        }
        int nMarkers = markers.getSize();
        for (int m = 0; m < nMarkers; m++)
        {
            wellMarker = (StsWellMarker) markers.getElement(m);
            StsModelSurface surface = wellMarker.getMarker().getModelSurface();
            if (velocityModel.hasModelSurface(surface)) continue;
            if (velocityModel.debug)
                oldTime = wellMarker.getLocation().getT();
            if (!computeMarkerTimeFromMDepth(wellMarker, velocityModel))
            {
                StsMessageFiles.errorMessage("Failed to convert well marker " + wellMarker.getName() + "; no time available.");
                continue;
            }
            if (velocityModel.debug)
                System.out.println("    well marker " + wellMarker.getName() + " well " + getName() + " readjusted from time " + oldTime + " to " + wellMarker.getLocation().getT());
        }
    }

    private void debugCheckMarkers()
    {
        StsObjectRefList markers = getMarkers();
        if (markers == null) return;
        int nMarkers = markers.getSize();
        for (int m = 0; m < nMarkers; m++)
        {
            StsWellMarker wellMarker = (StsWellMarker) markers.getElement(m);
            StsPoint location = wellMarker.getLocation();
            float markerT = location.getT();
            StsPoint markerPointOnWell = getPointAtZ(location.getZ(), false);
            float wellT = markerPointOnWell.getT();
            float timeError = markerT - wellT;
            //            if(timeError > 1.0f || timeError < -1.0f)
            System.out.println("DEBUG. Marker on surface for well " + getName() + " marker " + wellMarker.getName() + " time error " + timeError +
                " marker " + markerT + " well " + wellT);
        }
    } // after velocity model has been updated, compute depth values from time values

    // used for planned wells which are fixed in time
    public boolean adjustDepthPoints(StsSeismicVelocityModel velocityModel)
    {
        StsPoint[] points = getRotatedPoints();
        if (points == null) return false;
        for (int n = 0; n < points.length; n++)
        {
            float z;
            try
            {
                float[] xyzmt = points[n].v;
                z = (float) velocityModel.getZ(xyzmt[0], xyzmt[1], xyzmt[4]);
            }
            catch (Exception e)
            {
                StsMessageFiles.errorMessage("Failed to adjust well depth points for well " + getName());
                return false;
            }
            points[n].setZ(z);
        }
        return true;
    }

    public boolean canExport() { return true; }

    public boolean export()
    {
        return export(StsParameters.TD_ALL_STRINGS[zDomainSupported]);
    }

    public boolean export(String timeOrDepth)
    {
        return StsWellExportDialog.exportWell(currentModel, currentModel.win3d, "Well Export Utility", true, this, timeOrDepth);
    }


    /*
        private void outputPoint(StsPoint point0, StsPoint point1, float mDepth, boolean writeTime, boolean writeDepth, boolean exportLogData, StsObjectRefList curves, StsAsciiFile asciiFile)
        {
            String valLine = null;

            float m0 = point0.getM();
            float m1 = point1.getM();
            float x0 = point0.getX();
            float x1 = point1.getX();
            float y0 = point0.getY();
            float y1 = point1.getY();
            float z0 = point0.getZ();
            float z1 = point1.getZ();
            float t0 = point0.getT();
            float t1 = point1.getT();
            float f = (mDepth - m0)/(m1 - m0);
            float x = x0 + f*(x1 - x0);
            float y = y0 + f*(y1 - y0);
            float z = z0 + f*(z1 - z0);
            float t = t0 + f*(t1 - t0);
            float m = m0 + f*(m1 - m0);
            if (writeTime && writeDepth)
                 valLine = new String(x + " " + y + " " + z + " " + m + " " + t);
             else if (writeDepth)
                 valLine = new String(x + " " + y + z + " " + m);
             else if (writeTime)
                 valLine = new String(x + " " + y + " " + t);

            if(exportLogData)
            {
                for(int i=0; i<getNLogCurves(); i++)
                    valLine = valLine + " " + ((StsLogCurve)curves.getElement(i)).interpolatedValue(z);
            }
            try
            {
                asciiFile.writeLine(valLine);
            }
            catch(Exception e)
            {
            }
        }

        private void outputPoint(StsPoint point, boolean writeTime, boolean writeDepth, boolean exportLogData, StsObjectRefList curves, StsAsciiFile asciiFile)
        {
            String valLine = null;

            float m = point.getM();
            float x = point.getX();
            float y = point.getY();
            float z = point.getZ();
            float t = point.getT();
            if (writeTime && writeDepth)
                 valLine = new String(x + " " + y + " " + z + " " + m + " " + t);
             else if (writeDepth)
                 valLine = new String(x + " " + y + z + " " + m);
             else if (writeTime)
                 valLine = new String(x + " " + y + " " + t);
             if(exportLogData)
             {
                 for(int i=0; i<getNLogCurves(); i++)
                     valLine = valLine + " " + ((StsLogCurve)curves.getElement(i)).interpolatedValue(z);
             }
             try
             {
                 asciiFile.writeLine(valLine);
             }
             catch(Exception e)
             {
             }
         }
    */
    /*
        public void setExportName(String exportName)
        {
            this.exportName = exportName;
        }

        public String getExportName()
        {
            return exportName;
        }
    */
    public float getWellDirectionAngle()
    {
        StsPoint topPoint = getTopLineVertex().getPoint();
        StsPoint botPoint = getBotLineVertex().getPoint();
        float deltaX = topPoint.getX() - botPoint.getX();
        float deltaY = topPoint.getY() - botPoint.getY();
        return StsMath.atan2(deltaY, deltaX);
    }

    public boolean hasVsp() { return hasVsp; }

    public StsSeismicCurtain getSeismicCurtain() { return seismicCurtain; }

    public StsSeismicCurtain getCreateSeismicCurtain()
    {
        if (seismicCurtain == null)
        {
            isDrawingCurtain = true;
            createSeismicCurtain();
            currentModel.win3dDisplay();
        }
        return seismicCurtain;
    }

    public float[] getXYZForVertex(int vertexNum)
    {
        StsPoint[] points = getRotatedPoints();
        if (points == null) return null;
        if (points.length < vertexNum) return null;
        return points[vertexNum].getXYZorT();
    }

    public StsSeismicCurtain getCreateSeismicCurtain(StsSeismicVolume vol)
    {
        if (seismicCurtain == null)
            createSeismicCurtain(vol);
        return seismicCurtain;
    }

    public void setBornDate(String born)
    {
        super.setBornDate(born);
        bornField.setValueFromPanelObject(this);
        currentModel.viewObjectRepaint(this, this);
        return;
    }

    public void setDeathDate(String death)
    {
        super.setDeathDate(death);
        deathField.setValueFromPanelObject(this);
        currentModel.viewObjectRepaint(this, this);
        return;
    }

    // StsMonitorable Interface
    public int addNewData(StsObject object)
    {
        StsMessageFiles.errorMessage("addNewObject is not supportted for StsWell type.");
        return 0;
    }

    public int addNewData(double[] attValues, long time, String[] attNames)
    {
        StsMessageFiles.errorMessage("addNewData is not supportted for StsWell type.");
        return 0;
    }

    public int addNewData(StsPoint point, long time, String[] attNames)
    {
        StsMessageFiles.errorMessage("addNewData is not supportted for StsWell type.");
        return 0;
    }

    public int addNewData(String source, byte sourceType, long lastPollTime, boolean compute, boolean reload, boolean replace)
    {
        StsMessageFiles.errorMessage("addNewData is not supportted for StsWell type.");
        return 0;
    }

    public StsPerforationMarker[] getPerforationMarkers()
    {
        StsPerforationMarker[] list = null;
        if (markers != null)
        {
            int nMarkers = markers.getSize();
            for (int n = 0; n < nMarkers; n++)
            {
                StsWellMarker marker = (StsWellMarker) markers.getElement(n);
                if (marker instanceof StsPerforationMarker)
                    list = (StsPerforationMarker[])StsMath.arrayAddElement(list, marker);
            }
        }
        return list;
    }

    /* current marker file output example
    1G-03
    CURVE
    DEPTH
    VALUE
    TUZB 2710.790000
    WS2 3020.890000
    TWSK 3255.270000
    TWZD 3255.270000
    */

    public XMLObject getXMLobject()
    {
        return new XMLObject();
    }

    class XMLObject
    {
        String wellname;
        String depthType;
        StsWellMarker.XMLobject[] markerObjects;

        XMLObject()
        {
            wellname = name;
            depthType = StsWellKeywordIO.MDEPTH;
            int nMarkers = markers.getSize();
            {
                markerObjects = new StsWellMarker.XMLobject[nMarkers];
                for (int n = 0; n < nMarkers; n++)
                {
                    StsWellMarker marker = (StsWellMarker) markers.getElement(n);
                    markerObjects[n] = marker.getXMLobject();
                }
            }
        }
    }

    public StsAlarm[] getAlarms() { return alarms; }
    public boolean addAlarm(StsAlarm alarm)
    {
        if(alarms == null)
        {
            alarms = new StsAlarm[1];
            alarms[0] = alarm;
        }
        else
            alarms = (StsAlarm[])StsMath.arrayAddElement(alarms, alarm);

        return true;
    }
    public void checkAlarms()
    {
        for(int i=0; i<alarms.length; i++)
        {
            ;
        }
    }
    public boolean hasAlarms()
    {
        if(alarms != null)
            return true;
        else
            return false;
    }


    public void showPopupMenu(StsGLPanel glPanel, StsMouse mouse)
    {
        JPopupMenu tp = new JPopupMenu(name + " Properties");
        glPanel.add(tp);

        JLabel title = new JLabel("Well Popup Menu - " + name, JLabel.CENTER);
        title.setFont(new java.awt.Font("Dialog", 1, 14));
        tp.add(title);
        tp.addSeparator();
        
        StsMenuItem objectListWindow = new StsMenuItem();
        objectListWindow.setMenuActionListener("Group List...", this, "popupClassListPanel", null);
        tp.add(objectListWindow);

        if(Main.viewerOnly)
        {
            StsMenuItem classPropertyWindow = new StsMenuItem();
            classPropertyWindow.setMenuActionListener("Group Properties...", this, "popupClassPropertyPanel", null);
            tp.add(classPropertyWindow);

            StsMenuItem propertyWindow = new StsMenuItem();
            propertyWindow.setMenuActionListener("Display Properties...", this, "popupPropertyPanel", null);
            tp.add(propertyWindow);
            tp.addSeparator();
        }

        StsMenuItem wellWindow = new StsMenuItem();
        wellWindow.setMenuActionListener("Track Window...", this, "createWellWindow", null);
        tp.add(wellWindow);
        // StsMenuItem curtainView = new StsMenuItem();
        boolean hasObjects = false;
        ArrayList list = currentModel.getComboBoxToolbarPlugInDescriptors();
        for(int n = 0; n < list.size(); n++)
        {
            StsWorkflowPlugIn.ComboBoxDescriptor descriptor = (StsWorkflowPlugIn.ComboBoxDescriptor)list.get(n);
            if(descriptor.classNames == null) continue;
            hasObjects = true;
        }
        if(hasObjects)
        {
            StsMenuItem curtainWindow = new StsMenuItem();
            curtainWindow.setMenuActionListener("Curtain Window...", this, "addCurtainWindow", (StsGLPanel3d)glPanel);
            tp.add(curtainWindow);
        }
        tp.show(glPanel, mouse.getX(), mouse.getY());
    }

    public void createWellWindow()
    {
        openOrPopWindow();
    }

    public void addCurtainWindow(StsGLPanel3d glPanel3d)
    {
        StsObject object = ((StsComboBoxToolbar)currentModel.win3d.getToolbarNamed(StsComboBoxToolbar.NAME)).getSelectedObject();
        if(object instanceof StsVirtualVolume)
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Seismic curtains only support seismic volume, not virtual volumes.\n\nVirtual volume can be exported and re-imported as a seismic volume.");
        else
        {
            StsSeismicCurtainView seismicCurtainView = new StsSeismicCurtainView(currentModel, this);
            if(seismicCurtainView == null) return;
            StsWin3dBase window = glPanel3d.window;
            currentModel.createAuxWindow(window, StsSeismicCurtainView.shortViewNameCurtain, seismicCurtainView);
        }


     }

    private void addCurtainView(StsGLPanel3d glPanel3d)
    {
    }

}
