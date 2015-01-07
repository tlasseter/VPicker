package com.Sts.UI.Beans;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */


public class StsIntFieldBean extends StsLongFieldBean
{
    public StsIntFieldBean()
    {
    }

    public StsIntFieldBean(boolean editable)
    {
        super(editable, null);
    }

    public StsIntFieldBean(boolean editable, String label)
    {
        super(editable, label);
    }

    public StsIntFieldBean(Class c, String fieldName)
    {
        super(c, fieldName);
    }

    public StsIntFieldBean(Class c, String fieldName, String label)
    {
        super(c, fieldName, label);
    }

    public StsIntFieldBean(Class c, String fieldName, boolean editable, String label)
    {
        super(c, fieldName, editable, label);
    }

    public StsIntFieldBean(Class c, String fieldName, boolean editable)
    {
        super(c, fieldName, editable);
    }

    public StsIntFieldBean(Class c, String fieldName, int min, int max)
    {
        super(c, fieldName, min, max);
    }

    public StsIntFieldBean(Class c, String fieldName, int min, int max, String label)
    {
        super(c, fieldName, min, max, label);
    }

	public StsIntFieldBean(Class c, String fieldName, int min, int max, String label, boolean useArrows)
	{
		super(c, fieldName, min, max, label, useArrows);
    }

	public StsIntFieldBean(Class c, String fieldName, boolean editable, String label, boolean useArrows)
	{
		super(c, fieldName, editable, label, useArrows);
	}

    public StsIntFieldBean(Object beanObject, String fieldName, int min, int max, String label)
    {
        super(beanObject, fieldName, min, max, label);
    }

	public StsIntFieldBean(Object beanObject, String fieldName, int min, int max, String label, boolean useArrows)
	{
		super(beanObject, fieldName, min, max, label, useArrows);
    }

    public StsIntFieldBean(Object beanObject, String fieldName, boolean editable)
    {
        super(beanObject, fieldName, editable);
    }

	public StsIntFieldBean(Object beanObject, String fieldName, boolean editable, String label)
	{
		super(beanObject, fieldName, editable, label, false);
	}

	public StsIntFieldBean(Object beanObject, String fieldName, boolean editable, String label, boolean useArrows)
	{
		super(beanObject, fieldName, editable, label, useArrows);
	}

	public StsIntFieldBean(Object beanObject, String fieldName, String label, boolean useArrows)
	{
		super(beanObject, fieldName, true, label, useArrows);
	}

    public StsIntFieldBean(String labelText, int value)
    {
        setLabel(labelText);
        labelAndComponentLayout(label, textField);
        setValue(value);
        setEditable(false);
    }

	public StsIntFieldBean copy(Object beanObject)
	{
		StsIntFieldBean beanCopy = new StsIntFieldBean();
		if(editable)
			beanCopy.initialize(beanObject, fieldName, minimum, maximum, this.getLabelText(), useArrows);
		else
			beanCopy.initialize(beanObject, fieldName, false, getLabelText());
		return beanCopy;
	}

    public String toString()
	{
		return Integer.toString((int)value);
	}

	public Object fromString(String string)
	{
		return Integer.valueOf(string);
	}

	public Object getValueObject()
	{
		return new Integer((int)value);
	}
/*
	public Object getDefaultValueObject()
	{
		return new Integer((int)0);
	}
*/
	public int getIntValue()
	{
		 return Integer.parseInt(textField.getText());
	}

}
