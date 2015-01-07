package com.Sts.UI;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */

public class StsFilterProperties extends StsPanelProperties implements StsSerializable
{
    transient public final static byte NONE = 0;
    transient public final static byte PRESTACK = 1;
    transient public final static byte POSTSTACK = 2;
    transient public final static byte ALL_SEISMIC = 3;
    transient public final static byte VELOCITY = 4;
    transient public final static byte EVERYTHING = 5;
    transient static String[] whereStrings = new String[] {"None", "PreStack3d", "PostStack", "Pre & PostStack", "Velocity", "Seismic & Velocity"};
    transient static String[] velStrings = new String[] {"None", "Velocity" };

    // Box
	int filterWindowWidth = 5;
    byte whereBox = NONE;

    // Butterworth
	int bwLow= 0;
	int bwHigh = 50;
	int order = 2;
    byte whereBW = NONE;

    // Convolution
    byte whereConvolve = NONE;
    transient StsComboBoxFieldBean convKernelCombo;
    // Rank
    byte whereRank = NONE;
    double filterRadius = 3;
    transient StsComboBoxFieldBean rankSubTypeCombo;

    transient StsComboBoxFieldBean whereBoxBean;
    transient StsComboBoxFieldBean whereBWBean;
    transient StsComboBoxFieldBean whereConvolveBean;
    transient StsComboBoxFieldBean whereRankBean;

//	boolean rescaleRequired = true; // non-transient so it can be propagated

	transient StsJPanel panel = null;
	transient StsGroupBox filterBox = null;

	transient StsGroupBox filtBox = null;

	static private final String title = "Filter Properties";
	public StsFilterProperties()
	{
	}

	public StsFilterProperties(String fieldName)
	{
		super(title, fieldName);
    }

	public StsFilterProperties(StsObject parentObject, StsFilterProperties defaultProperties, String fieldName)
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
        whereBoxBean = new StsComboBoxFieldBean(this, "whereBoxString", "Apply Box To:", whereStrings);
        whereBWBean = new StsComboBoxFieldBean(this, "whereBWString", "Apply Butterworth To:", whereStrings);

        whereConvolveBean = new StsComboBoxFieldBean(this, "whereConvolveString", "Apply Convolution To:", velStrings);
//        convKernelCombo = new StsComboBoxFieldBean(this,"kernelString","Kernel:", StsConvolve.KERNELS);

        whereRankBean = new StsComboBoxFieldBean(this, "whereRankString", "Apply Rank To:", velStrings);
//        rankSubTypeCombo = new StsComboBoxFieldBean(this,"subTypeString","Sub-Type:", StsRankFilters.RANKFILTERS);

		propertyBeans = new StsFieldBean[]
		{
            whereBoxBean,
			new StsIntFieldBean(this, "filterWindowWidth", 3, 33, "Box Filter Width", true),
            whereBWBean,
            new StsIntFieldBean(this, "bwLow", 0, 250, "Lower Frequency (Hz)", true),
			new StsIntFieldBean(this, "bwHigh", 0, 250, "Upper Frequency (Hz)", true),
			new StsIntFieldBean(this, "order", 2, 8, "Filter Order (2-8)", true),
            whereConvolveBean,
//            convKernelCombo,
            whereRankBean,
//            rankSubTypeCombo,
            new StsDoubleFieldBean(this,"filterRadius", 1, 20, "Radius:", true)
		};
	}

	public int getFilterWindowWidth() { return filterWindowWidth; }

	public boolean getApplyBoxFilter(byte dataType)
    {
        return isApplied(dataType, whereBox);
    }
	public boolean getApplyBWFilter(byte dataType)
    {
        return isApplied(dataType, whereBW);
    }
    public boolean getApplyConvolveFilter(byte dataType)
    {
        return isApplied(dataType, whereConvolve);
    }
    public boolean getApplyRankFilter(byte dataType)
    {
        return isApplied(dataType, whereRank);
    }

    public boolean isApplied(byte dataType, byte whereApplied)
    {
        if(whereApplied == NONE) return false;
        switch(dataType)
        {
            case PRESTACK:
                if((whereApplied == PRESTACK) || (whereApplied == ALL_SEISMIC) || (whereApplied == EVERYTHING))
                    return true;
                else
                    return false;
            case POSTSTACK:
                if((whereApplied == POSTSTACK) || (whereApplied == ALL_SEISMIC) || (whereApplied == EVERYTHING))
                    return true;
                else
                    return false;
            case VELOCITY:
                if((whereApplied == VELOCITY) || (whereApplied == EVERYTHING))
                    return true;
                else
                    return false;
            default:
                break;
        }
        return false;
    }

    public byte getWhereBoxApplied() { return whereBox; }
    public byte getWhereBWApplied() { return whereBW; }
    public byte getWhereConvolveApplied() { return whereConvolve;  }
    public byte getWhereRankApplied()
    {
        return whereRank;
    }

    public String getWhereBoxString() { return whereStrings[whereBox]; }
    public String getWhereBWString() { return whereStrings[whereBW]; }

    // Index is from Where list even though this is not the list used in UI for Convolve and Rank
    public String getWhereConvolveString()
    {
        return whereStrings[whereConvolve];
    }
    public String getWhereRankString()
    {
        return whereStrings[whereRank];
    }

    public void setWhereBoxString(String where)
    { 
        whereBox = whereApplied(where);
    }
    public void setWhereBWString(String where) 
    { 
        whereBW = whereApplied(where);
    }
    public void setWhereConvolveString(String where) 
    { 
        whereConvolve = velApplied(where);
    }
    public void setWhereRankString(String where)
    { 
        whereRank = velApplied(where); 
    }

    public byte whereApplied(String stg)
    {
        for(int i=0; i<whereStrings.length; i++)
        {
            if (stg.equals(whereStrings[i]))
                return (byte)i;
        }
        return NONE;
    }
    public byte velApplied(String stg)
    {
        byte where = NONE;
        for(int i=0; i<velStrings.length; i++)
        {
            if (stg.equals(velStrings[i]))
                where = (byte)i;
        }
        if(where == 1)
            return VELOCITY;
        return NONE;
    }
    // Buterworth Filter Properties
	public int getOrder() { return order; }
	public int getBwLow() { return bwLow; }
	public int getBWHigh() { return bwHigh; }

    // Box Filter Properties
	public void setFilterWindowWidth(int width)
	{ 
	    filterWindowWidth = 1+ 2*(width/2);
	}
	public void setBwLow(int val) 
	{ 
	    bwLow = val;
	}
	public void setBwHigh(int val)
	{ 
	    bwHigh = val;
	}
	public void setOrder(int val)
	{ 
	    order = val;
	}

    // Rank Filter Properties
    public void setSubTypeString(String subType) { }
 //   public String getSubTypeString() { return StsRankFilters.RANKFILTERS[getSubType()]; }
    public byte getSubType()
    {
        return (byte)rankSubTypeCombo.getSelectedIndex();
    }
    public double getFilterRadius() { return filterRadius; }
    public void setFilterRadius(double radius)
    {
        filterRadius = radius;
    }
    // Convolution Filter Properties
    public void setKernelString(String kernel) { }
//    public String getKernelString() { StsConvolve.KERNELS[getKernel()];}
    public byte getKernel()
    {
        return (byte)convKernelCombo.getSelectedIndex();
    }
	static int getStringIndex(String string, String[] strings)
	{
		for(int n = 0; n < strings.length; n++)
			if(strings[n] == string)return n;
		return -1;
	}

	static public int getColorIndex(Color color, Color[] colors)
	{
		for(int n = 0; n < colors.length; n++)
			if(color.equals(colors[n]))return n;
		return 0;
	}
}
