package com.Sts.UI.ObjectPanel;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class ObjectPanelTest implements StsTreeObjectI
{
    protected boolean booleanCheck1 = false;
	protected boolean booleanCheck2 = false;
    protected float floatValue = 100.0f;
    protected int intValue = -10;
    StsColorscale colorscale = new StsColorscale("Test", StsSpectrum.createRainbowSpectrum(StsSpectrumClass.SPECTRUM_RAINBOW, 255), -100.0f, 100.0f);
//    protected Color color = Color.white;

    static StsFieldBean[] displayFields = null;
    static StsFieldBean[] propertyFields = null;

    public ObjectPanelTest()
    {
    }

    public boolean getBooleanCheck1() { return booleanCheck1; }
    public void setBooleanCheck1(boolean check) {booleanCheck1 = check; }

	public boolean getBooleanCheck2() { return booleanCheck2; }
	public void setBooleanCheck2(boolean check) {booleanCheck2 = check; }

    public float getFloatBox() { return floatValue; }
    public void setFloatBox(float value) { floatValue = value; }

    public int getIntBox() { return intValue; }
    public void setIntBox(int value) { intValue = value; }

//    public Color getColor() { return color; }
//    public void setBeachballColors(Color color) {stsColor = color; }


	public void setColorscale(StsColorscale colorscale)
	{
		this.colorscale = colorscale;
	}

	public StsColorscale getColorscale()
	{
		return colorscale;
	}

    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields == null)
        {
            displayFields = new StsFieldBean[]
            {
                new StsBooleanFieldBean(ObjectPanelTest.class, "booleanCheck1", "Boolean Check 1"),
                new StsBooleanFieldBean(ObjectPanelTest.class, "booleanCheck2", "Boolean Check 2"),
                new StsFloatFieldBean(ObjectPanelTest.class, "floatBox", "Float Box"),
                new StsIntFieldBean(ObjectPanelTest.class, "intBox", "Int Box")
            };
        }
        return displayFields;
    }

    public StsFieldBean[] getDefaultFields() { return null; }
    public StsFieldBean[] getPropertyFields()
    {
        if(propertyFields == null)
        {
            propertyFields = new StsFieldBean[]
            {
                new StsEditableColorscaleFieldBean(ObjectPanelTest.class, "colorscale")
            };
        }
        return propertyFields;
    }

    public Object[] getChildren() { return new Object[0]; }
    public StsObjectPanel getObjectPanel() { return null; }
    public String getName() { return "Test Object"; }
    public void treeObjectSelected() { }

    public boolean export() { return false; }
    public boolean launch() { return false; }
    public boolean canExport() { return false; }
    public boolean canLaunch() { return false; }
    public void popupPropertyPanel() { return; }
   
//    public boolean makeCurrent() { return false; }
    public boolean anyDependencies() { return false; }

    void printState()
    {
        System.out.println("CheckBox1: " + booleanCheck1);
		System.out.println("CheckBox2: " + booleanCheck2);
        System.out.println("FloatValue: " + floatValue);
        colorscale.debugPrint();
    }
}
