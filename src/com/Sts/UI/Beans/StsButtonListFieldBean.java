
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Beans;

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2005</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author TJLasseter
 * @version c62e
 *
 * A series of JRadioButtons each with an associated buttonObject are ganged
 * together in a ButtonGroup, i.e., only one can be selected at a time.  Each
 * buttonObject must implement a readable toString() method; this string is the
 * label for the button.
 * When a button is selected, the associated beanObject to which this bean is
 * connected has it's setter called with the field name of this bean.  See and
 * run main() method for an example of usage.
 */
public class StsButtonListFieldBean extends StsFieldBean
{
	/** prefix which contains all buttons added, one button for each object in the items list */
	ButtonGroup group = new ButtonGroup();
	/** set of buttons: one button for each buttonObject */
    JRadioButton[] buttons = null;
	/** buttonObjects: one object for each button */
	Object[] buttonObjects = null;
	/** Panel which contains this button prefix: can be added as a single component to a panel */
	StsJPanel buttonPanel = StsJPanel.addInsets();

	public StsButtonListFieldBean()
    {
    }

	public StsButtonListFieldBean(Class c, String fieldName, String fieldLabel, Object[] items, boolean isHorizontal)
    {
        initialize(c, fieldName, fieldLabel, items, isHorizontal);
    }

    public StsButtonListFieldBean(Object instance, String fieldName, String fieldLabel, Object[] items, boolean isHorizontal)
    {
        initialize(instance, fieldName, fieldLabel, items, isHorizontal);
    }

    public void initialize(Class c, String fieldName, String fieldLabel, Object[] items, boolean isHorizontal)
    {
		setListItems(items);
		classInitialize(c, fieldName, fieldLabel);
//		setLabel(fieldLabel);
		layoutBean(isHorizontal);
        addActionListener();
    }

    public void initialize(Object instance, String fieldName, String fieldLabel, Object[] items, boolean isHorizontal)
    {
		setListItems(items);
        super.initialize(instance, fieldName, true, fieldLabel);
        layoutBean(isHorizontal);
        addActionListener();
    }

	public StsButtonListFieldBean copy(Object beanObject)
	{
		StsButtonListFieldBean beanCopy = new StsButtonListFieldBean();
		beanCopy.initialize(beanObject, fieldName, getLabelText(), buttonObjects, false);
		return beanCopy;
	}

	protected void layoutBean(boolean isHorizontal)
	 {
		 try
		 {
			 setLayout(new GridBagLayout());
			 GridBagConstraints gbc = new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0);

			 if (label != null) add(label, gbc);
			 gbc.gridx++;

			 layoutButtons(isHorizontal);
			 add(buttonPanel, gbc);

			 validate();
		 }
		 catch (Exception e)
		 {
			 e.printStackTrace();
		 }
	}

	private void layoutButtons(boolean isHorizontal)
	{
		buttonPanel.gbc.anchor = GridBagConstraints.WEST;
		if(isHorizontal)
		{
			if(label != null) buttonPanel.addToRow(label);
			for(int n = 0; n < buttons.length; n++)
				buttonPanel.addToRow(buttons[n]);
		}
		else
		{
			if(label != null) buttonPanel.add(label);
			for (int n = 0; n < buttons.length; n++)
				buttonPanel.add(buttons[n]);
		}
	}

    public void addActionListener()
    {
		for( int i=0; i<buttons.length; i++ )
        {
        	buttons[i].addActionListener(this);
//            buttons[i].setActionCommand(getLabel());
        }
    }

    public Object getValueObject()
	{
		for(int i=0; i<buttons.length; i++ )
			if(buttons[i].isSelected()) return buttonObjects[i];
		return null;
    }
/*
	public Object getDefaultValueObject()
	{
		if(buttonObjects != null && buttonObjects.length > 0)
			return buttonObjects[0];
		else
			return null;
	}
*/
    public String toString()
    {
		for(int i=0; i<buttons.length; i++ )
        	if( buttons[i].isSelected() )
            	return buttons[i].getText();
        return null;
    }

	public Object fromString(String string)
	{
		for(int i=0; i<buttons.length; i++ )
			if(buttons[i].getText().equals(string))
            	return buttons[i];
		return null;
	}

    public void setListItems(Object[] items)
    {
    	if( items == null || items.length == 0) return;
		buttons = new JRadioButton[items.length];
		buttonObjects = items;
		for( int i=0; i<items.length; i++ )
		{
			buttons[i] = new JRadioButton(items[i].toString());
			group.add(buttons[i]);
        }
//		buttons[0].setSelected(true);
	}

    public void doSetValueObject(Object object)
    {
        if(buttons == null) return;
		for(int i=0; i<buttonObjects.length; i++ )
        {
        	if(buttonObjects[i] == object)
            	buttons[i].setSelected(true);
        }
    }

    public Component[] getBeanComponents() { return new Component[] { label, buttonPanel }; }
	public Component getButtonPanel() { return buttonPanel; }

    public void actionPerformed(ActionEvent e)
    {
        if(buttons == null || buttons.length == 0) return;
        setValueInPanelObjects();
    }


	protected void setEditable()
	{
		buttonPanel.setEditable(editable);
	}

	static public class TestButtonList
	{
		String buttonASelected = buttonASelections[1];
		TestButtonSelection buttonBSelected = buttonBSelections[0];
		boolean checkBox1 = false;
		boolean checkBox2 = true;

	    static String[] buttonASelections = new String[] { "Button-A-One", "Button-A-Two" };
		static TestButtonSelection buttonB1 = new TestButtonSelection("Button-B-One");
		static TestButtonSelection buttonB2 = new TestButtonSelection("Button-B-Two");

	    static TestButtonSelection[] buttonBSelections = new TestButtonSelection[] { buttonB1, buttonB2 };

		TestButtonList()
		{
		}

		public String getButtonASelected() { return buttonASelected; }
		public void setButtonASelected(String buttonSelected) { this.buttonASelected = buttonSelected; }
		public TestButtonSelection getButtonBSelected() { return buttonBSelected; }
		public void setButtonBSelected(TestButtonSelection buttonSelected) { this.buttonBSelected = buttonSelected; }
		public void setCheckBox1(boolean value) { checkBox1 = value; }
		public boolean getCheckBox1() { return checkBox1; }
		public void setCheckBox2(boolean value) { checkBox2 = value; }
		public boolean getCheckBox2() { return checkBox2; }
	}

	static public class TestButtonSelection
	{
		String name;
		TestButtonSelection(String name)
		{
			this.name = name;
		}

		public String toString() { return name; }
	}

	static public void main(String[] args)
	{
		JFrame frame = new JFrame("RadioButton List Test");
		 frame.setSize(300, 200);

		 Container contentPane = frame.getContentPane();
		 StsGroupBox groupBox = new StsGroupBox("Select-O-Button");
		 TestButtonList testButtonList = new TestButtonList();
		 StsButtonListFieldBean buttonAList = new StsButtonListFieldBean(testButtonList, "buttonASelected", null, testButtonList.buttonASelections, true);
		 StsButtonListFieldBean buttonBList = new StsButtonListFieldBean(testButtonList, "buttonBSelected", null, testButtonList.buttonBSelections, true);
		 groupBox.gbc.gridwidth = 2;
		 groupBox.add(buttonAList);
		 groupBox.add(buttonBList);
		 StsBooleanFieldBean checkBoxBean1 = new StsBooleanFieldBean(testButtonList, "checkBox1", "Checkbox 1");
		 StsBooleanFieldBean checkBoxBean2 = new StsBooleanFieldBean(testButtonList, "checkBox2", "Checkbox 2");
		 groupBox.gbc.gridwidth = 1;
		 groupBox.addToRow(checkBoxBean1);
		 groupBox.addEndRow(checkBoxBean2);
	     contentPane.add(groupBox);
		 StsToolkit.centerComponentOnScreen(frame);
		 frame.pack();
         frame.setVisible(true);
	}
}
