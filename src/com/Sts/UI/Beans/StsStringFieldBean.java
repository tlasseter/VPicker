package com.Sts.UI.Beans;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

public class StsStringFieldBean extends StsTextFieldBean
{
    String value = "";

    public StsStringFieldBean()
    {
    }

    public StsStringFieldBean(String label)
    {
        this(true, label);
    }

    public StsStringFieldBean(boolean editable)
    {
        this(editable, null);
    }

    public StsStringFieldBean(boolean editable, String label)
    {
        this.classInitialize((Class)null, null, editable, label);
    }

    public StsStringFieldBean(String fieldName, String label)
    {
        this.initialize(null, fieldName, value, true, label);
    }

    public StsStringFieldBean(Class c, String fieldName)
    {
        this.initialize(c, fieldName, value, true, null);
    }

    public StsStringFieldBean(Class c, String fieldName, String label)
    {
        this.initialize(c, fieldName, value, true, label);
    }

    public StsStringFieldBean(Class c, String fieldName, boolean editable, String label)
    {
        this.initialize(c, fieldName, value, editable, label);
    }

    public StsStringFieldBean(Class c, String fieldName, boolean editable)
    {
        this.initialize(c, fieldName, value, editable, "");
    }

    public StsStringFieldBean(Class c, String fieldName, String value, String label)
    {
        this.initialize(c, fieldName, value, true, label);
    }

	public StsStringFieldBean(Class c, String fieldName, String value, boolean editable, String label)
	{
		this.initialize(c, fieldName, value, editable, label);
	}

	public StsStringFieldBean(String value, boolean editable, String label)
	{
		this.initialize(null, null, value, editable, label);
	}

	public StsStringFieldBean(Object beanObject, String fieldName, boolean editable)
	{
		this(beanObject, fieldName, editable, null);
		layoutBean();
    }

    public StsStringFieldBean(Object beanObject, String fieldName, String label)
	{
		this(beanObject, fieldName, true, label);
		layoutBean();
    }

    public StsStringFieldBean(Object beanObject, String fieldName, boolean editable, String label)
    {
        initialize(beanObject, fieldName, editable, label);
		layoutBean();
    }

    public StsStringFieldBean(Object beanObject, String fieldName, String value, boolean editable, String label)
    {
        this.initialize(beanObject, fieldName, value, editable, label);
    }
/*
    public void classInitialize(String label)
    {
        classInitialize(null, null, "", true, label);
    }
*/
    public void initialize(Class c, String fieldName, String value, String label)
    {
        super.classInitialize(c, fieldName, true, label);
		initializeValue(value);
		layoutBean();
   }

   private void initializeValue(String value)
   {
	   this.value = value;
	   super.setText(value);
   }

    public void initialize(Class c, String fieldName, String value, boolean editable, String label)
    {
        super.classInitialize(c, fieldName, editable, label);
        initializeValue(value);
		layoutBean();
    }

	public void classInitialize(Class c, String fieldName, boolean editable, String label)
	{
		super.classInitialize(c, fieldName, editable, label);
		layoutBean();
	}

    public void initialize(Object beanObject, String fieldName, String value, boolean editable, String label)
    {
        super.initialize(beanObject, fieldName, editable, label);
        initializeValue(value);
		layoutBean();
    }

	public void initialize(Object beanObject, String fieldName, boolean editable, String label)
	{
		super.initialize(beanObject, fieldName, editable, label);
		layoutBean();
    }

	public StsStringFieldBean copy(Object beanObject)
	{
		StsStringFieldBean beanCopy = new StsStringFieldBean();
		beanCopy.initialize(beanObject, fieldName, editable, getLabelText());
		return beanCopy;
	}

	protected void layoutBean()
	{
		textField.setColumns(10);
		labelAndComponentLayout(label, textField);
    }

    public void setTextForeground(Color color)
    {
        textField.setForeground(color);
    }

    public void setTextFont(Font font)
    {
        textField.setFont(font);
    }

    public JTextField getTextField()
    {
        return textField;
    }
    public void setTextBackground(Color color)
    {
        textField.setBackground(color);
    }
/*
	protected void layoutBean()
	{
		try
		{
			if (label != null)
			{
				gbc.weightx = 0.0;
				gbc.fill = GridBagConstraints.NONE;
				addToRow(label);
				gbc.fill = GridBagConstraints.HORIZONTAL;
				gbc.weightx = 1.0;
			}
			addToRow(textField);
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
    }
*/
    // return true if value has changed
    public boolean setValue(String newValue)
    {
        try
        {
            String oldValue = value;
            value = newValue;
            setText(newValue);
            return !newValue.equals(oldValue);
        }
        catch(Exception e)
        {
            StsToolkit.beep();
            textField.requestFocus();
            return false;
        }
    }

	public String toString() { return value; }
	public Object fromString(String string) { return string; }

    public void doSetValueObject(Object object)
    {
        if(!(object instanceof String)) return;
        setValue((String)object);
    }

    public Object getValueObject()
    {
        return value;
    }

//	public Object getDefaultValueObject() { return new String(""); }

    public String getValue()
    {
        return textField.getText();
    }

    protected boolean isStringBean() { return true; }
}
