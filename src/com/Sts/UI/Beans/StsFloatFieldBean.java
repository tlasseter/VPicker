package com.Sts.UI.Beans;

/**
 * <p>Title: Field Beans Development</p>
 * <p>Description: General beans for generic panels.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.0
 */

import com.Sts.Utilities.*;

import java.awt.*;

public class StsFloatFieldBean extends StsDoubleFieldBean
{
    public StsFloatFieldBean()
	{
    }

    public StsFloatFieldBean(boolean editable)
    {
        super(editable, null);
    }

    public StsFloatFieldBean(boolean editable, String label)
    {
        super(editable, label);
    }

    public StsFloatFieldBean(Class c, String fieldName)
    {
        super(c, fieldName);
    }

    public StsFloatFieldBean(Class c, String fieldName, String label)
    {
        super(c, fieldName, label);
    }

    public StsFloatFieldBean(Class c, String fieldName, boolean editable, String label)
    {
        super(c, fieldName, editable, label);
    }

    public StsFloatFieldBean(Class c, String fieldName, boolean editable)
    {
        super(c, fieldName, editable);
    }

    public StsFloatFieldBean(Class c, String fieldName, float min, float max)
    {
        super(c, fieldName, min, max);
    }

    public StsFloatFieldBean(Class c, String fieldName, float min, float max, String label)
    {
        super(c, fieldName, min, max, label);
    }

	public StsFloatFieldBean(Class c, String fieldName, float min, float max, String label, boolean useArrows)
	{
		super(c, fieldName, true, label, useArrows);
		setRange(min, max);
    }

    public StsFloatFieldBean(Class c, String fieldName, boolean editable, String label, boolean useArrows)
	{
		super(c, fieldName, editable, label, useArrows);
    }

    public StsFloatFieldBean(Object beanObject, String fieldName, String label)
    {
        super(beanObject, fieldName, label);
    }

    public StsFloatFieldBean(Object beanObject, String fieldName)
    {
        super(beanObject, fieldName);
    }

    public StsFloatFieldBean(Object beanObject, String fieldName, boolean editable)
    {
        super(beanObject, fieldName, editable, null);
    }

	public StsFloatFieldBean(Object beanObject, String fieldName, boolean editable, String label)
	{
		super(beanObject, fieldName, editable, label);
	}

	public StsFloatFieldBean(Object beanObject, String fieldName, boolean editable, String label, boolean useArrows)
	{
		super(beanObject, fieldName, editable, label, useArrows);
    }

	public StsFloatFieldBean(Object beanObject, String fieldName, float min, float max, String label)
	{
		super(beanObject, fieldName, min, max, label);
    }

	public StsFloatFieldBean(Object beanObject, String fieldName, float min, float max, String label, boolean useArrows)
	{
		super(beanObject, fieldName, min, max, label, useArrows);
    }

    public StsFloatFieldBean(String labelText, float value)
    {
        setLabel(labelText);
        textField.setText(formatNumber(value));
        labelAndComponentLayout(label, textField);
        setEditable(false);
    }

	public StsFloatFieldBean copy(Object beanObject)
	{
		StsFloatFieldBean beanCopy = new StsFloatFieldBean();
		if(editable)
			beanCopy.initialize(beanObject, fieldName, minimum, maximum, this.getLabelText(), useArrows);
		else
			beanCopy.initialize(beanObject, fieldName, false, getLabelText());
		return beanCopy;
	}

    public Object getValueObject()
	{
		return new Float(value);
    }

	public String toString()
	{
		return Float.toString((float)value);
	}
	public Object fromString(String string)
	{
		return Float.valueOf(string);
	}
/*
	public Object getDefaultValueObject()
	{
		return new Float(0);
	}
*/
	public float getFloatValue()
	{
		return (float)value;
	}

	public String formatNumber(double number)
	{
		return formatNumber(number, 6);
	}

	public boolean isValidString(String string)
	{
		if(string == null || string.length() <= 0) return true;
		try
		{
			float newValue = labelFormat.parse(string).floatValue();
//			double newValue = labelFormat.pDouble.parseDouble(string);
//			setText(string);
//			value = newValue;
			return true;
		}
		catch(Exception e)
		{
			setValue(value);
			return false;
		}
    }

    static public void main(String[] args)
    {
		TestFloatBean testFloatBean = new TestFloatBean();
		StsFloatFieldBean floatBean1 = new StsFloatFieldBean(testFloatBean, "value1", true, "Value 1", true);
		floatBean1.setValueAndRangeFixStep(1.0, 0.0f, 1000000.0f, 1000.0f);
		StsJPanel panel = new StsJPanel();
		panel.add(floatBean1);
        panel.setMaximumSize(new Dimension(100, 100));
        StsToolkit.createDialog(panel, true, 100, 100);
        floatBean1.resetLabel("Changed");
		testFloatBean.print();
		floatBean1.print();
    }
}
/** Note that we can have explicit get/set (as for value1) or not (value2).  */
class TestFloatBean
{
    private float value1 = 1000.0f;
	private float value2 = 10.0f;
	private float value3 = 123456.789f;
    public void setValue1(float v) { value1 = v; System.out.println("Value set to " + v); }
    public float getValue1() { return value1; }
	public void setValue3(float v) { value3 = v; System.out.println("Value set to " + v);  }
    public float getValue3() { return value3; }
	public void print() { System.out.println("Bean panel values. value1: " + value1 + " value2: " + value2 + " value3: " + value3); }
}
