package com.Sts.UI.Beans;

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: Field Beans Development</p>
 * <p>Description: General beans for generic panels.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.0
 */

public class StsRadioButtonFieldBean extends StsFieldBean
{
    String selectedText;
    String unselectedText = null;
    JRadioButton radioButton = new JRadioButton(selectedText);

    public StsRadioButtonFieldBean()
    {
    }

    public StsRadioButtonFieldBean(Object beanObject, String fieldName, String fieldLabel, ButtonGroup group)
    {
        initialize(beanObject, fieldName, fieldLabel, group);
    }
    
    public void initialize(Object beanObject, String fieldName, String fieldLabel)
    {
        this.beanObject = beanObject;
        Class c = beanObject.getClass();
        selectedText = fieldLabel;
        radioButton.setText(fieldLabel);
        classInitialize(c, fieldName, fieldLabel);
        setValueFromPanelObject(beanObject);
    }

    public void initialize(Object beanObject, String fieldName, String fieldLabel, ButtonGroup group)
	{
        checkLabel(fieldLabel);
        selectedText = fieldLabel;
        this.beanObject = beanObject;
		Class c = beanObject.getClass();
		selectedText = fieldLabel;
		radioButton.setText(fieldLabel);
        layoutBean();
        classInitialize(c, fieldName, fieldLabel);
		setValueFromPanelObject(beanObject);
        group.add(radioButton);
        addActionListener();
    }

	public StsRadioButtonFieldBean copy(Object beanObject)
	{
		StsRadioButtonFieldBean beanCopy = new StsRadioButtonFieldBean();
		beanCopy.initialize(beanObject, fieldName, getLabelText());
		return beanCopy;
	}

    private void checkLabel(String fieldLabel)
    {
        if(fieldLabel == null || fieldLabel.length() == 0)
            StsException.systemError(this, "checkLabel", "checkBox cannot be null or zero-length string");
    }

    protected void layoutBean()
    {
        try
        {
			add(radioButton);
        }
        catch(Exception e) { e.printStackTrace(); }
    }

    public void setSelected(boolean isSelected)
    {
        radioButton.setSelected(isSelected);
        if(unselectedText == null) return;
        if(isSelected) radioButton.setText(selectedText);
        else           radioButton.setText(unselectedText);
        repaint();
    }

	public String toString()
	{
		return label + " " + Boolean.toString(radioButton.isSelected());
	}

	public Object fromString(String string)
	{
		return Boolean.valueOf(string);
	}

    public boolean isSelected() { return radioButton.isSelected(); }
    public void setText(String text)
    {
        selectedText = text;
        radioButton.setText(text);
    }
    public String getText() { return selectedText; }
    public void setUnselectedText(String text) { unselectedText = text; }
    public String getUnselectedText() { return unselectedText; }

    public void setHorizontalAlignment(int alignment) { radioButton.setHorizontalAlignment(alignment); }
    public void setHorizontalTextPosition(int alignment) { radioButton.setHorizontalTextPosition(alignment); }

    public void addActionListener()
    {
        radioButton.addActionListener(this);
        radioButton.setActionCommand(fieldName);
    }

    public void doSetValueObject(Object valueObject)
    {
        if(!(valueObject instanceof Boolean)) return;
        boolean selected = ((Boolean)valueObject).booleanValue();
        radioButton.setSelected(selected);
    }

    public Object getValueObject()
    {
        return new Boolean(radioButton.isSelected());
    }
/*
	public Object getDefaultValueObject()
	{
		return new Boolean(true);
	}
*/
	public void actionPerformed(ActionEvent e)
	{
		setValueInPanelObjects();
	}

	public void setValueInPanelObjects()
	{
		Object[] panelObjects = getPanelObjects();
		if (panelObjects == null)return;
		Object panelObject = null;
		Object value = selectedText;

		try
		{
			if (set != null)
			{
				for (int n = 0; n < panelObjects.length; n++)
				{
					panelObject = panelObjects[n];
					if (panelObject != null) set.invoke(panelObject, value);
				}
			}
			else if (field != null)
			{
				for (int n = 0; n < panelObjects.length; n++)
				{
					panelObject = panelObjects[n];
					if (panelObject != null) field.set(panelObject, value);
				}
			}
		}
		catch (Exception e)
		{
			outputFieldException(e, panelObject, value);
		}
	}

    public void addButtonToGroup(ButtonGroup group)
    {
        group.add(radioButton);
    }

    public Component[] getBeanComponents() { return new Component[] { label, radioButton }; }
	public Component getMainComponent() { return radioButton; }

    public static void main(String[] args)
    {
        StsRadioButtonFieldBean stsRadioButtonFieldBean1 = new StsRadioButtonFieldBean();
    }
}