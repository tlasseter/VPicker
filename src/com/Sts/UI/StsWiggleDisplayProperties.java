package com.Sts.UI;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.util.*;

/**
 * <p>Title: S2S development</p>
 * <p/>
 * <p>Description: Integrated seismic to simulation software</p>
 * <p/>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p/>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */

public class StsWiggleDisplayProperties extends StsSeismicPanelProperties implements StsSerializable
{
    protected StsColor wigglePlusColor = CLEAR;
    protected StsColor wiggleMinusColor = CLEAR;
    protected StsColor wiggleBackgroundColor = StsColor.WHITE;
    protected boolean wiggleDrawLine = true;
    protected boolean wiggleDrawPoints = false;
    protected boolean wiggleSmoothCurve = false;
    protected boolean wiggleReversePolarity = false;
    protected boolean drawAttributeCurve = false;

    protected boolean drawNMOCurves = true;
    
    protected int offsetAxisTypeIndex = OFFSET_AXIS_ABS_INDEX;
    protected int attributeIndex = -1; //, dAttributeIndex = 0, vAttributeIndex = 0;
    protected String attributeName = ATTRIBUTE_NONE; //, dAttributeIndex = 0, vAttributeIndex = 0;
    public float stretchMute = 50.0f;
    private float stretchMuteMinOffset = initStretchMuteMinOffset(); //even after initialization, this keeps getting set to zero ???
    public float inchesPerSecond = 3.f;
    public float tracesPerInch = 15.f;
    public enum DisplayTypes { WiggleTrace, VariableArea };
    public DisplayTypes displayType = DisplayTypes.WiggleTrace;
    {
        wiggleOverlapPercent = 30;
    }
    
//	transient float currentTracesPerInch = tracesPerInch;
//    transient float currentInchesPerSecond = inchesPerSecond;

    protected boolean displayWiggles = true;

//	public transient boolean axisTypeChanged = false;
//    public transient boolean axisRangeChanged = false;

    /** checked by Gather view on a viewObjectChanged(wiggleDisplayProperties) to see if axisType has changed */
    transient public boolean axisTypeChanged = false;
    transient public boolean stretchMuteChanged = false;
    transient StsModel model;

    transient String[] attributes = null;
//    transient String[] dAttributes = null;
//    transient String[] vAttributes = null;
//    transient StsComboBoxFieldBean datumSelectBean = null;
    //    transient StsComboBoxFieldBean velocitySelectBean = null;
    transient StsComboBoxFieldBean attributeSelectBean = null;
//	transient StsButtonFieldBean displayColorscaleButtonBean;
    //	transient StsButtonFieldBean editColorscaleButtonBean;
    transient StsFloatFieldBean tpiBean = null;
    transient StsFloatFieldBean ipsBean = null;
    transient StsIntFieldBean wiggleOverlapPercentBean = null;
    transient StsFloatFieldBean stretchMuteBean = null;
    transient StsFloatFieldBean stretchMuteMinOffsetBean = null;

    static public final String title = "Wiggle Display Properties";

	static public String ATTRIBUTE_NONE = "None";
    /** List of offset axis types */
    static public final int OFFSET_AXIS_VALUE = 0;
    static public final int OFFSET_AXIS_ABS_INDEX = 1;
    static public final int OFFSET_AXIS_INDEX = 2;
    static public final String OFFSET_AXIS_VALUE_STRING = "Value";
    static public final String OFFSET_AXIS_ABS_INDEX_STRING = "Absolute value index";
    static public final String OFFSET_AXIS_INDEX_STRING = "Value index";
    static public final String[] offsetAxisTypeStrings = new String[]
            {OFFSET_AXIS_VALUE_STRING, OFFSET_AXIS_ABS_INDEX_STRING, OFFSET_AXIS_INDEX_STRING};

    static public final StsColor CLEAR = new StsColor(0, 0, 0, 0);
    /** List of possible plus colors */
    static public final StsColor[] plusColors = new StsColor[]
            {StsColor.BLACK, CLEAR, StsColor.WHITE, StsColor.GRAY, StsColor.RED, StsColor.BLUE};
    /** List of possible minus colors */
    static public final StsColor[] minusColors = new StsColor[]
            {CLEAR, StsColor.BLACK, StsColor.WHITE, StsColor.GRAY, StsColor.RED, StsColor.BLUE};
    /** List of possible background colors; line will be other color */
    static public final StsColor[] backgroundColors = new StsColor[]
            {StsColor.WHITE, StsColor.BLACK, StsColor.GRAY, StsColor.LIGHT_GRAY, StsColor.DARK_GRAY};
    static public final StsColor[] lineColors = new StsColor[]
            {StsColor.BLACK, StsColor.WHITE, StsColor.WHITE, StsColor.BLACK, StsColor.WHITE};

    static final public float stretchMuteMax = 1000.0f;
    static final public int overlapPercentMax = 10000; //sometimes higher gain needed SWC 6/8/09

    public StsWiggleDisplayProperties()
    {
        //initializeBeans();  jbw defer building until after isVSP is restored.
//        System.out.println("Default WiggleDisplayProperties constructor");
    }

    private float initStretchMuteMinOffset() {
    	StsPreStackLineSet lineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
    	if (lineSet != null) {
    		StsPreStackLine[] lines = lineSet.lines;
    		if (lines != null && lines.length > 0) {
    			StsPreStackLine line = lines[0];
    			if (line != null) {
    				stretchMuteMinOffset = line.cdpAvgInterval * 4;
    			}
    		}
    	}
    	return stretchMuteMinOffset;
    }

	public StsWiggleDisplayProperties(StsClass stsClass, String fieldName)
    {
        super(title, fieldName);
        parentObject = stsClass;
//        this.isVSP = isVSP;
    }

	public StsWiggleDisplayProperties(StsClass stsClass, String fieldName, String attributeName)
    {
        super(title, fieldName);
        parentObject = stsClass;
		this.attributeName = attributeName;
//        this.isVSP = isVSP;
    }

    public StsWiggleDisplayProperties(Object parentObject, StsWiggleDisplayProperties defaultProperties, String fieldName)
    {
        super(parentObject, defaultProperties, title, fieldName);
        initializeDefaultProperties(defaultProperties);
//        this.isVSP = isVSP;
    }

    public void initializeDefaultProperties(Object defaultProperties)
    {
        if(defaultProperties != null)
           StsToolkit.copyObjectNonTransientFields(defaultProperties, this);
    }
       /*
      public void saveState()
      {
       savedProperties = (StsWiggleDisplayProperties)StsToolkit.copyObjectNonTransients(this);
      }
      */
    public void buildBeanLists()
    {
        attributes = new String[]{ATTRIBUTE_NONE};
        if (parentObject == null) return;

        String[] volAttributes = null;
        if(parentObject instanceof StsPreStackLineSetClass)
        {

            StsPreStackLineSet volume = StsPreStackLineSetClass.currentProjectPreStackLineSet;
            /*
            if (volume.lines == null) return;

            if (volume.getZDomain() == StsParameters.TD_TIME)
                volAttributes = volume.lines[0].getTimeAttributes();
            else
                volAttributes = volume.lines[0].getDistanceAttributes();
                */
            volAttributes = volume.attributeNames; //why not let user have full range of selection? Besides, OFFSET is most important and wasn't in list SWC 12/3/09
            
        }
        else if (parentObject instanceof StsVsp)
        {
            volAttributes = ((StsVsp) parentObject).getTimeAttributes();
            volAttributes = (String[]) StsMath.arrayAddArray(volAttributes, ((StsVsp) parentObject).getDistanceAttributes());
        }
        if (volAttributes != null)
        {
            Arrays.sort(volAttributes);
            attributes = (String[]) StsMath.arrayAddArray(attributes, volAttributes);
        }
    }

    public void initializeBeans()
    {
        buildBeanLists();

        attributeSelectBean = new StsComboBoxFieldBean(this, "AttributeName", "Display Trace Attribute:", attributes);
        tpiBean = new StsFloatFieldBean(this, "tracesPerInch", 1.f, 250.f, "Traces per Inch:");
        ipsBean = new StsFloatFieldBean(this, "inchesPerSecond", .1f, 50.0f, "Inches per Second:");
        wiggleOverlapPercentBean = new StsIntFieldBean(this, "wiggleOverlapPercent", 0, overlapPercentMax, "Gain (Overlap Percent):", true);
        wiggleOverlapPercentBean.setStep(25);
        stretchMuteBean = new StsFloatFieldBean(this, "stretchMute", 0, stretchMuteMax, "Stretch Mute (%):");
        stretchMuteBean.fixStep(1.0);
        stretchMuteMinOffsetBean = new StsFloatFieldBean(this, "stretchMuteMinOffset", 0, 100000, "Min Mute Offset:", true);
        stretchMuteMinOffsetBean.fixStep(1.0);
        if (!wiggleDrawLine) displayType = DisplayTypes.VariableArea;
        if (parentObject instanceof StsPreStackLineSetClass)
        {
//			velocitySelectBean = new StsComboBoxFieldBean(this, "velocityAttribute", "Velocity for Datum Correction:", vAttributes);
//			datumSelectBean = new StsComboBoxFieldBean(this, "datumAttribute", "Correct Traces to Datum:", dAttributes);
            propertyBeans = new StsFieldBean[]
            {
//			    new StsBooleanFieldBean(this, "displayWiggles", "Display wiggles"),
                    tpiBean,
                    ipsBean,
                    new StsColorComboBoxFieldBean(this, "wigglePlusColor", "Plus Color", plusColors),
                    new StsColorComboBoxFieldBean(this, "wiggleMinusColor", "Minus Color", minusColors),
                    new StsColorComboBoxFieldBean(this, "wiggleBackgroundColor", "Background Color", backgroundColors),
                    wiggleOverlapPercentBean,
                    stretchMuteBean,
                    stretchMuteMinOffsetBean,
                    new StsComboBoxFieldBean(this, "displayType", "Display type:", DisplayTypes.values()),
                    //new StsBooleanFieldBean(this, "wiggleDrawLine", "Variable Area"),
                    new StsBooleanFieldBean(this, "wiggleReversePolarity", "Reverse Polarity"),
                    new StsBooleanFieldBean(this, "wiggleDrawPoints", "Draw Data Points"),
                    new StsBooleanFieldBean(this, "wiggleSmoothCurve", "Smooth Curve    "),
                    new StsBooleanFieldBean(this, "drawNMOCurves", "Draw NMO Curves    "),
                    new StsComboBoxFieldBean(this, "offsetAxisTypeString", "Offset Axis Type", offsetAxisTypeStrings),
                    attributeSelectBean,
                    new StsBooleanFieldBean(this, "drawAttributeCurve", "Plot Attribute Over Traces"),
            };
        }
        else if(parentObject instanceof StsVsp)
        {
            propertyBeans = new StsFieldBean[]
            {
                new StsBooleanFieldBean(this, "displayWiggles", "Display wiggles"),
                tpiBean,
                ipsBean,
                new StsColorComboBoxFieldBean(this, "wigglePlusColor", "Plus Color", plusColors),
                new StsColorComboBoxFieldBean(this, "wiggleMinusColor", "Minus Color", minusColors),
                new StsColorComboBoxFieldBean(this, "wiggleBackgroundColor", "Background Color", backgroundColors),
                wiggleOverlapPercentBean,
                new StsBooleanFieldBean(this, "wiggleDrawLine", "Draw Trace Line "),
                new StsBooleanFieldBean(this, "wiggleReversePolarity", "Reverse Polarity"),
                attributeSelectBean,
            };
        }
        if (parentObject instanceof StsSeismicVolumeClass)
        {
            propertyBeans = new StsFieldBean[]
            {
                    new StsColorComboBoxFieldBean(this, "wigglePlusColor", "Plus Color", plusColors),
                    new StsColorComboBoxFieldBean(this, "wiggleMinusColor", "Minus Color", minusColors),
                    new StsColorComboBoxFieldBean(this, "wiggleBackgroundColor", "Background Color", backgroundColors),
                    wiggleOverlapPercentBean,
                    new StsBooleanFieldBean(this, "wiggleDrawLine", "Draw Trace Line  "),
                    new StsBooleanFieldBean(this, "wiggleDrawPoints", "Draw Data Points"),
                    new StsBooleanFieldBean(this, "wiggleSmoothCurve", "Smooth Curve    ")
            };
        }
    }

    public void setWigglePlusColor(StsColor color)
    {
        this.wigglePlusColor = color;
    }

    public void setWiggleMinusColor(StsColor color)
    {
        wiggleMinusColor = color;
    }

    public void setWiggleBackgroundColor(StsColor color)
    {
        wiggleBackgroundColor = color;
    }

    public void setWiggleDrawLine(boolean draw)
    {
        wiggleDrawLine = draw;
    }

    public void setWiggleDrawPoints(boolean draw)
    {
        wiggleDrawPoints = draw;
    }

    public void setWiggleSmoothCurve(boolean smoothCurve)
    {
        wiggleSmoothCurve = smoothCurve;
    }

    public void setWiggleReversePolarity(boolean reversePolarity)
    {
        wiggleReversePolarity = reversePolarity;
    }

    public float getTracesPerInch() { return tracesPerInch; }

    public float getInchesPerSecond() { return inchesPerSecond; }

    public void setTracesPerInch(float tracesPerInch)
    {
        if (this.tracesPerInch == tracesPerInch) return;
        this.tracesPerInch = tracesPerInch;
        setRangeChanged(true);
        if (tpiBean != null) tpiBean.setValue(tracesPerInch);
    }

    public void setInchesPerSecond(float inchesPerSecond)
    {
        if (this.inchesPerSecond == inchesPerSecond) return;
        this.inchesPerSecond = inchesPerSecond;
        setRangeChanged(true);
        if (ipsBean != null) ipsBean.setValue(inchesPerSecond);
    }

    public void setOffsetAxisTypeString(String type)
    {
        int typeIndex = getStringIndex(type, offsetAxisTypeStrings);
        if (typeIndex == offsetAxisTypeIndex) return;
        offsetAxisTypeIndex = typeIndex;
        axisTypeChanged = true;
    }

    public void setAttributeName(String type)
    {
        attributeName = type;
        //attributeSelectBean.setSelectedItem(attributeName);
        
        // Set bean to value if required
        if(attributes == null)
            return;
        attributeIndex = getStringIndex(type, attributes);
        /*
        if (attributeSelectBean.getSelectedIndex() != index)
            attributeSelectBean.setSelectedIndex(index);
        */
            
    }

    /**
     * Turn wiggle display on and off. When on it will use the wiggle to pixel ratio to determine when
     * wigles are displayed. Currently only works for 2D displays
     *
     * @param displayWiggles boolean
     */
    public void setDisplayWiggles(boolean displayWiggles)
    {
        if (this.displayWiggles == displayWiggles) return;
        this.displayWiggles = displayWiggles;
    }

    public StsColor getWigglePlusColor()
    {
        return this.wigglePlusColor;
    }

    public StsColor getWiggleMinusColor()
    {
        return wiggleMinusColor;
    }

    public StsColor getWiggleBackgroundColor()
    {
        return wiggleBackgroundColor;
    }

    public boolean getWiggleDrawLine()
    {
        return wiggleDrawLine;
    }

    public boolean getWiggleDrawPoints()
    {
        return wiggleDrawPoints;
    }

    public boolean getWiggleSmoothCurve()
    {
        return wiggleSmoothCurve;
    }

    public boolean getWiggleReversePolarity()
    {
        return wiggleReversePolarity;
    }

    public boolean getDisplayWiggles()
    {
        return displayWiggles;
    }

    public StsColor getLineColor()
    {
        for (int n = 0; n < backgroundColors.length; n++)
            if (wiggleBackgroundColor.equals(backgroundColors[n]))
            {
                return lineColors[n];
            }
        return StsColor.WHITE;
    }

    public int getOffsetAxisType()
    {
        return offsetAxisTypeIndex;
    }

    public String getOffsetAxisTypeString()
    {
        return offsetAxisTypeStrings[offsetAxisTypeIndex];
    }

    public String getAttributeName()
    {
        return attributeName;
        /*
        if (attributes == null)
            return "None";
        else
            return attributes[attributeIndex];
         */
    }

    public float getStretchMute() { return stretchMute; }

    public void setStretchMute(float mute)
    {
        stretchMute = mute;
        stretchMuteChanged = true;
        if (this.stretchMuteBean != null) stretchMuteBean.setValue(stretchMute);
    }
    
    public float getStretchMuteMinOffset() { return stretchMuteMinOffset; }

    public void setStretchMuteMinOffset(float offset)
    {
    	stretchMuteMinOffset = offset;
        stretchMuteChanged = true;
        if (this.stretchMuteMinOffsetBean != null) stretchMuteMinOffsetBean.setValue(offset);
    }

    /*
      public void setTracesRescaled()
      {
       isChanged = false;
      }

      public boolean getSetTracesRescaleRequired()
      {
       if(!isChanged) return false;
       isChanged = false;
       return true;
      }
      */
    public boolean hasFill()
    {
        // return false;
        return wigglePlusColor != CLEAR || wiggleMinusColor != CLEAR;
    }

    public void displayWiggleProperties()
    {
        StsOkCancelDialog dialog = new StsOkCancelDialog(currentModel.win3d, this, "Edit Wiggle Display Properties", false);
    }

    public void displayWiggleProperties(String title)
    {
        StsOkCancelDialog dialog = new StsOkCancelDialog(currentModel.win3d, this, title, false);
    }

    public void dialogSelectionType(int selectionType)
    {
        super.dialogSelectionType(selectionType);
        axisTypeChanged = false;
        stretchMuteChanged = false;
    }

    static int getStringIndex(String string, String[] strings)
    {
        for (int n = 0; n < strings.length; n++)
            if (strings[n] == string) return n;
        return -1;
    }

    static public int getColorIndex(Color color, Color[] colors)
    {
        for (int n = 0; n < colors.length; n++)
            if (color.equals(colors[n])) return n;
        return 0;
    }
    
    public void setDisplayType( DisplayTypes type )
    {
        if (type.equals(DisplayTypes.WiggleTrace))
        {
            wiggleDrawLine = true;
            displayType = DisplayTypes.WiggleTrace;
        }
        else
        {
            wiggleDrawLine = false;
            displayType = DisplayTypes.VariableArea;
        }
    }
/*
	public void displayColorscale()
	{
		StsColorscalePanel colorPanel = new StsColorscalePanel(((StsVsp)parentObject).getColorscale(), true);
		StsToolkit.createDialog(colorPanel);
	}

	public void editColorscale()
	{
		StsColorscalePanel colorPanel = new StsColorscalePanel(((StsVsp)parentObject).getColorscale(), true);
		currentModel.getActionManager().startAction(StsColorscaleAction.class, new Object[] { colorPanel, (Frame)currentModel.win3d } );
//		redraw = true;
		currentModel.win3dDisplayAll();
	}
*/
    
    public boolean getDrawNMOCurves()
    {
        return drawNMOCurves;
    }

    public void setDrawNMOCurves(boolean drawNMOCurves)
    {
        this.drawNMOCurves = drawNMOCurves;
    }


    public boolean getDrawAttributeCurve()
    {
        return drawAttributeCurve;
    }

    public void setDrawAttributeCurve(boolean drawAttributeCurve)
    {
        this.drawAttributeCurve = drawAttributeCurve;
    }
}