package com.Sts.UI;

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

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;


public class StsMarkerProperties extends StsPanelProperties
{
    protected int[] colorIdx = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
    transient StsColorComboBoxFieldBean colorListBean[] = new StsColorComboBoxFieldBean[32];

	transient public StsModel model;

	static private final String title = "Marker Properties";

	public StsMarkerProperties()
	{
	}

	public StsMarkerProperties(StsModel model, StsClass stsClass, String fieldName)
	{
		super(title, fieldName);
		this.model = model;
	}

	public StsMarkerProperties(StsObject parentObject, StsMarkerProperties defaultProperties, String fieldName)
	{
        super(parentObject, title, fieldName);
        initializeDefaultProperties(defaultProperties);
    }
    public void initializeDefaultProperties(Object defaultProperties)
    {
        if(defaultProperties != null)
            StsToolkit.copyObjectNonTransientFields(defaultProperties, this);
    }

	public void buildBeanLists(Object parentObject)
	{
        for(int i=0; i<colorIdx.length; i++)
        {
            colorListBean[i] = new StsColorComboBoxFieldBean(this, "color"+i, "Marker #" + i + ": ", StsColor.colors32);
            //colorListBean[i].setBeachballColors(StsColor.stsColors32[colorIdx[i]]);
        }
    }

    public void initializeBeans()
    {
        buildBeanLists(parentObject);
        propertyBeans = colorListBean;
    }
    public StsColor getColor(int index)
    {
        index = index%colorIdx.length;
        return StsColor.colors32[colorIdx[index]];
    }

    public StsColor getColor0() { return StsColor.colors32[colorIdx[0]]; }
    public StsColor getColor1() { return StsColor.colors32[colorIdx[1]]; }
    public StsColor getColor2() { return StsColor.colors32[colorIdx[2]]; }
    public StsColor getColor3() { return StsColor.colors32[colorIdx[3]]; }
    public StsColor getColor4() { return StsColor.colors32[colorIdx[4]]; }
    public StsColor getColor5() { return StsColor.colors32[colorIdx[5]]; }
    public StsColor getColor6() { return StsColor.colors32[colorIdx[6]]; }
    public StsColor getColor7() { return StsColor.colors32[colorIdx[7]]; }
    public StsColor getColor8() { return StsColor.colors32[colorIdx[8]]; }
    public StsColor getColor9() { return StsColor.colors32[colorIdx[9]]; }
    public StsColor getColor10() { return StsColor.colors32[colorIdx[10]]; }
    public StsColor getColor11() { return StsColor.colors32[colorIdx[11]]; }
    public StsColor getColor12() { return StsColor.colors32[colorIdx[12]]; }
    public StsColor getColor13() { return StsColor.colors32[colorIdx[13]]; }
    public StsColor getColor14() { return StsColor.colors32[colorIdx[14]]; }
    public StsColor getColor15() { return StsColor.colors32[colorIdx[15]]; }
    public StsColor getColor16() { return StsColor.colors32[colorIdx[16]]; }
    public StsColor getColor17() { return StsColor.colors32[colorIdx[17]]; }
    public StsColor getColor18() { return StsColor.colors32[colorIdx[18]]; }
    public StsColor getColor19() { return StsColor.colors32[colorIdx[19]]; }
    public StsColor getColor20() { return StsColor.colors32[colorIdx[20]]; }
    public StsColor getColor21() { return StsColor.colors32[colorIdx[21]]; }
    public StsColor getColor22() { return StsColor.colors32[colorIdx[22]]; }
    public StsColor getColor23() { return StsColor.colors32[colorIdx[23]]; }
    public StsColor getColor24() { return StsColor.colors32[colorIdx[24]]; }
    public StsColor getColor25() { return StsColor.colors32[colorIdx[25]]; }
    public StsColor getColor26() { return StsColor.colors32[colorIdx[26]]; }
    public StsColor getColor27() { return StsColor.colors32[colorIdx[27]]; }
    public StsColor getColor28() { return StsColor.colors32[colorIdx[28]]; }
    public StsColor getColor29() { return StsColor.colors32[colorIdx[29]]; }
    public StsColor getColor30() { return StsColor.colors32[colorIdx[30]]; }
    public StsColor getColor31() { return StsColor.colors32[colorIdx[31]]; }

     public void setColor0(StsColor color) { colorIdx[0] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor1(StsColor color) { colorIdx[1] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor2(StsColor color) { colorIdx[2] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor3(StsColor color) { colorIdx[3] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor4(StsColor color) { colorIdx[4] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor5(StsColor color) { colorIdx[5] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor6(StsColor color) { colorIdx[6] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor7(StsColor color) { colorIdx[7] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor8(StsColor color) { colorIdx[8] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor9(StsColor color) { colorIdx[9] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor10(StsColor color) { colorIdx[10] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor11(StsColor color) { colorIdx[11] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor12(StsColor color) { colorIdx[12] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor13(StsColor color) { colorIdx[13] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor14(StsColor color) { colorIdx[14] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor15(StsColor color) { colorIdx[15] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor16(StsColor color) { colorIdx[16] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor17(StsColor color) { colorIdx[17] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor18(StsColor color) { colorIdx[18] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor19(StsColor color) { colorIdx[19] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor20(StsColor color) { colorIdx[20] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor21(StsColor color) { colorIdx[21] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor22(StsColor color) { colorIdx[22] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor23(StsColor color) { colorIdx[23] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor24(StsColor color) { colorIdx[24] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor25(StsColor color) { colorIdx[25] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor26(StsColor color) { colorIdx[26] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor27(StsColor color) { colorIdx[27] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor28(StsColor color) { colorIdx[28] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor29(StsColor color) { colorIdx[29] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor30(StsColor color) { colorIdx[30] = StsColor.getColorIndex(color, StsColor.colors32); }
     public void setColor31(StsColor color) { colorIdx[31] = StsColor.getColorIndex(color, StsColor.colors32); }
}