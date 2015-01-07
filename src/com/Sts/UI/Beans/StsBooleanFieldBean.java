package com.Sts.UI.Beans;

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

public class StsBooleanFieldBean extends StsFieldBean
{
	JCheckBox checkBox = new JCheckBox();
	private boolean labelOnLeft = false;
//	StsJPanel panel = StsJPanel.noInsets();

	public StsBooleanFieldBean()
	{
		super();
		layoutBean();
		addActionListener();
	}

	public StsBooleanFieldBean(Class c, String fieldName, String label)
	{
		super();
		classInitialize(c, fieldName, label);
		addActionListener();
	}

	public StsBooleanFieldBean(Object beanObject, String fieldName, String label)
	{
		super();
		initialize(beanObject, fieldName, label);
		addActionListener();
	}

	public StsBooleanFieldBean(Object beanObject, String fieldName, String label, boolean labelOnLeft)
	{
		super();
		initialize(beanObject, fieldName, label, labelOnLeft);
		addActionListener();
	}

	public StsBooleanFieldBean(Class c, String fieldName, boolean initialValue, String label)
	{
		super();
		classInitialize(c, fieldName, label);
		checkBox.setSelected(initialValue);
		addActionListener();
	}

	public StsBooleanFieldBean(Class c, String fieldName, String label, boolean editable, boolean labelOnLeft)
	{
		super();
        this.labelOnLeft = labelOnLeft;
        classInitialize(c, fieldName, editable, label);
        layoutBean();
        addActionListener();
	}

    public StsBooleanFieldBean(Class c, String fieldName, String label, boolean editable)
	{
		super();
		classInitialize(c, fieldName, editable, label);
        layoutBean();
        addActionListener();
	}

    public StsBooleanFieldBean(String label)
	{
		super();
		classInitialize(null, null, true, label);
        layoutBean();
	}

    public void initialize(Object beanObject, String fieldName, String fieldLabel)
	{
		super.initialize(beanObject, fieldName, true, fieldLabel);
		layoutBean();
		//addActionListener();
	}

	public void initialize(Object beanObject, String fieldName, String fieldLabel, boolean labelOnLeft)
	{
		super.initialize(beanObject, fieldName, true, fieldLabel);
		this.labelOnLeft = labelOnLeft;
		layoutBean();
		//addActionListener();
	}

	public void classInitialize(Class c, String fieldName, String fieldLabel)
	{
		super.classInitialize(c, fieldName, fieldLabel);
		layoutBean();
		//addActionListener();
	}

	public StsBooleanFieldBean(Object beanObject, String fieldName, boolean initialValue, String label)
	{
		super();
		initialize(beanObject, fieldName, label);
		checkBox.setSelected(initialValue);
		addActionListener();
	}

	public StsBooleanFieldBean(Object beanObject, String fieldName, boolean initialValue, String label, boolean labelOnLeft)
	{
		super();
		this.labelOnLeft = labelOnLeft;
		initialize(beanObject, fieldName, label);
		checkBox.setSelected(initialValue);
		addActionListener();
	}

	public StsBooleanFieldBean copy(Object beanObject)
	{
		StsBooleanFieldBean beanCopy = new StsBooleanFieldBean();
		beanCopy.initialize(beanObject, fieldName, getLabelText());
		return beanCopy;
	}
/*
	protected void layoutBean()
	{
		add(checkBox, BorderLayout.WEST);
		if(label != null)
			add(label, BorderLayout.EAST);
	}
*/
	protected void layoutBean()
	{
		try
		{
			checkBox.setMargin(new Insets(0, 0, 0, 2));
			gbc.anchor = GridBagConstraints.WEST;
//			gbc.fill = GridBagConstraints.NONE;

			if (labelOnLeft)
			{
				if (label != null)
				{
					gbc.weightx = 1.0;
					label.setHorizontalAlignment(SwingConstants.LEFT);
					label.setHorizontalTextPosition(SwingConstants.LEFT);
					addToRow(label);
				}
				gbc.weightx = 0.0;
				addToRow(checkBox);
			}
			else
			{
				gbc.weightx = 0.0;
				addToRow(checkBox);
				if (label == null)
					return;
				gbc.weightx = 1.0;
				label.setHorizontalAlignment(SwingConstants.LEFT);
				label.setHorizontalTextPosition(SwingConstants.LEFT);
				addToRow(label);
			}
			gbc.anchor = GridBagConstraints.CENTER;
//			gbc.fill = GridBagConstraints.HORIZONTAL;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	// no way that I can figure out how to make JCheckBox uneditable
	public void setEditable()
	{
		checkBox.setEnabled(editable);
//		checkBox.getModel().setEditable(editable);
	}

	/** Treat a check box and its label as a single component on any panel */
	public void addComponentsToPanel(StsJPanel panel)
	{
		panel.add(this, panel.gbc);
		panel.gbc.gridx += gbc.gridwidth;
	}

	public void setSelected(boolean isSelected) { checkBox.setSelected(isSelected); }
	public boolean isSelected() { return checkBox.isSelected(); }
	public void setHorizontalAlignment(int alignment) { checkBox.setHorizontalAlignment(alignment); }
	public void setHorizontalTextPosition(int alignment) { checkBox.setHorizontalTextPosition(alignment); }

	public void addActionListener()
	{
		checkBox.addActionListener(this);
		checkBox.setActionCommand(fieldName);
	}

	public String toString()
	{
		boolean isSelected = checkBox.isSelected();
		return Boolean.toString(isSelected);
	}

	public Object fromString(String string)
	{
		return Boolean.valueOf(string);
	}

	public void doSetValueObject(Object valueObject)
	{
		boolean selected;
		if(valueObject instanceof Boolean)
		{
			selected = ((Boolean)valueObject).booleanValue();
			checkBox.setSelected(selected);
			return;
		}
		else if(valueObject instanceof String)
		{
			if(valueObject.equals("true"))
				checkBox.setSelected(true);
			else
				checkBox.setSelected(false);
		}
	}

	public void setValue(boolean selected)
	{
        if(checkBox.isSelected() == selected) return;
        checkBox.setSelected(selected);
	}

	public Object getValueObject()
	{
		return new Boolean(checkBox.isSelected());
	}
/*
	public Object getDefaultValueObject()
	{
		return new Boolean(true);
	}
*/
	public boolean getValue() { return checkBox.isSelected(); }

	public Component[] getBeanComponents() { return new Component[] { this }; }
    public Component getMainComponent() { return checkBox; }
    public JCheckBox getCheckBox() { return checkBox; }
//    public Component[] getBeanComponents() { return new Component[] { checkBox, label }; }
	public void actionPerformed(ActionEvent e)
	{
		setValueInPanelObjects();
	}

	public Object decode(String string)
	{
		return new Boolean(string);
	}

	protected void copy()
	{
		checkBox = new JCheckBox();
		removeAll();
		layoutBean();
		addActionListener();
	}

    public static void main(String[] args)
	{
		StsBooleanFieldBean StsBooleanFieldBean1 = new StsBooleanFieldBean();
	}
}
