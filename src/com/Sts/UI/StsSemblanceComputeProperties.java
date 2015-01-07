package com.Sts.UI;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

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
public class StsSemblanceComputeProperties extends StsPanelProperties implements StsSerializable
{
    public byte semblanceType = SEMBLANCE_STANDARD;
    public int windowWidth = 15;
    public byte order = ORDER_FAST;
    public boolean normalizeCoherence = true;
    public double clip = 0.6;

    public boolean enhancedAdjust = false;

    transient public StsModel model;
    transient protected StsPreStackLineSet3dClass seismicClass;

    transient StsBooleanFieldBean enhancedAdjustBean;
    transient StsBooleanFieldBean normalizeCoherenceBean;
    transient StsDoubleFieldBean clipCoherenceBean;

    static public final String SEMBLANCE_STANDARD_STRING = "Semblance";
    static public final String SEMBLANCE_ENHANCED_STRING = "Enhanced Semblance";
    static public final String SEMBLANCE_STACKED_AMP_STRING = "Stacked Amplitude";
    //	static public final String SEMBLANCE_JIGGLE_STRING = "Enhanced plus jiggle";
    static final String[] semblanceTypes = new String[]{SEMBLANCE_STANDARD_STRING, SEMBLANCE_ENHANCED_STRING, SEMBLANCE_STACKED_AMP_STRING};
    static public final byte SEMBLANCE_STANDARD = 0;
    static public final byte SEMBLANCE_ENHANCED = 1;
    public static final byte SEMBLANCE_STACKED_AMP = 2;
//	static public final int SEMBLANCE_JIGGLE = 2;

    /** order of NMO equation for semblanceBytes and gathers */
    static public final String ORDER_2ND_STRING = "Normal (2nd)";
    static public final String ORDER_4TH_STRING = "4th Order";
    static public final String ORDER_6TH_STRING = "6th Order";
    static public final String ORDER_FAST_STRING = "Fast (2nd)";
    static public final String ORDER_6TH_OPT_STRING = "Opt 6th Order";
    static public final String[] orderStrings = new String[]{ORDER_2ND_STRING, ORDER_4TH_STRING, ORDER_6TH_STRING, ORDER_6TH_OPT_STRING, ORDER_FAST_STRING};
    static public final byte ORDER_2ND = 0;
    static public final byte ORDER_4TH = 1;
    static public final byte ORDER_6TH = 2;
    static public final byte ORDER_6TH_OPT = 3;
    static public final byte ORDER_FAST = 4;

    static public final String title = "Coherence Compute Properties";
    static public final boolean immediateChange = true;
    static public final String ClipCoherenceText = "Clip Colorscale: "; //don't know why, but space after colorscale: keeps motion of components after label change to a minimum
    private static final String MaxSemblanceText = "Max Semblance:";

    transient public boolean redraw = false;

    public StsSemblanceComputeProperties()
    {
    }

    public StsSemblanceComputeProperties(String fieldName)
    {
        super(title, fieldName);
    }

    public StsSemblanceComputeProperties(StsObject parentObject, StsSemblanceComputeProperties defaultProperties, String fieldName)
    {
        super(parentObject, title, fieldName);
        initializeDefaultProperties(defaultProperties);
    }

    public void initializeDefaultProperties(Object defaultProperties)
    {
        if(defaultProperties != null)
            StsToolkit.copyObjectNonTransientFields(defaultProperties, this);
    }

    public void initializeBeans()
    {
        model = StsObject.getCurrentModel();
        seismicClass = (StsPreStackLineSet3dClass)model.getStsClass(StsPreStackLineSet3d.class);
        enhancedAdjustBean = new StsBooleanFieldBean(this, "enhancedAdjust", "Enhanced adjust");
        enhancedAdjustBean.getLabel().setEnabled(false);
        enhancedAdjustBean.getCheckBox().setEnabled(false);
        normalizeCoherenceBean = new StsBooleanFieldBean(this, "normalizeCoherence", "Saturate Colorscale");
        clipCoherenceBean = new StsDoubleFieldBean(this, "clip", 0, 1, ClipCoherenceText, true);

        propertyBeans = new StsFieldBean[]
        {
                new StsButtonListFieldBean(this, "semblanceTypeString", "Coherence Types:", semblanceTypes, false),
                normalizeCoherenceBean,
                clipCoherenceBean,
                new StsIntFieldBean(this, "windowWidth", 0, 1000, "Window Width:", true),
                new StsComboBoxFieldBean(this, "orderString", "Semblance Order", orderStrings),
                enhancedAdjustBean
        };
    }

    public String getOrderString() { return orderStrings[order]; }

    public byte getOrder() { return order; }

    public String getSemblanceTypeString() { return semblanceTypes[semblanceType]; }

    public int getWindowWidth() { return windowWidth; }

    public void setOrder(byte order) 
    {
        this.order = order;
    }

    public void setSemblanceTypeString(String string)
    {
        semblanceType = getSemblanceTypeFromString(string);
        if(semblanceType == SEMBLANCE_ENHANCED)
        {
            enhancedAdjustBean.getLabel().setEnabled(true);
            enhancedAdjustBean.getCheckBox().setEnabled(true);
        }
        else
        {
            enhancedAdjustBean.getLabel().setEnabled(false);
            enhancedAdjustBean.getCheckBox().setEnabled(false);
        }
        if(semblanceType == SEMBLANCE_STANDARD || semblanceType == SEMBLANCE_ENHANCED)
        {
            normalizeCoherenceBean.getLabel().setEnabled(true);
            normalizeCoherenceBean.getCheckBox().setEnabled(true);
        }
        if(semblanceType == SEMBLANCE_STACKED_AMP)
        {
            normalizeCoherenceBean.getLabel().setEnabled(false);
            normalizeCoherenceBean.getCheckBox().setSelected(true);
            normalizeCoherenceBean.getCheckBox().setEnabled(false);
            setNormalizeCoherence(true);
        }
    }

    public void setWindowWidth(int width)
    { 
        windowWidth = width;
    }

    public void setOrderString(String orderString)
    {
        for(int n = 0; n < orderStrings.length; n++)
        {
            if(orderString == orderStrings[n])
            {
                this.order = (byte)n;
                return;
            }
        }
    }

    private byte getSemblanceTypeFromString(String typeString)
    {
        for(byte n = 0; n < semblanceTypes.length; n++)
            if(typeString == semblanceTypes[n]) return n;
        return 0;
    }

    public boolean isNormalizeCoherence()
    {
        return normalizeCoherence;
    }

    public void setNormalizeCoherence(boolean normalizeCoherence)
    {
        this.normalizeCoherence = normalizeCoherence;
        if(normalizeCoherence)
        {
            //clipCoherenceBean.resetLabel(ClipCoherenceText);
            clipCoherenceBean.getLabel().setText(ClipCoherenceText);
        }
        else
        {
            clipCoherenceBean.getLabel().setText(MaxSemblanceText);
        }
    }

    public boolean getEnhancedAdjust() { return enhancedAdjust; }

    public void setEnhancedAdjust(boolean adjust)
    {
        this.enhancedAdjust = adjust;
    }

    public double getClip()
    {
        return clip;
    }

    public void setClip(double clip)
    {
        this.clip = clip;
    }
}
