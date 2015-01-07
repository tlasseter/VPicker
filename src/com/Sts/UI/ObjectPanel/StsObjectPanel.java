
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.ObjectPanel;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class StsObjectPanel extends StsJPanel implements StsBeanPanelI
{
	StsFieldBean[] displayFields = null;
	StsFieldBean[] defaultFields = null;
    StsTreeObjectI[] panelObjects = null;
	StsFieldBean[] propertyFields = null;
	int[] selected = null;
//	Method getName = null;
//	Method setName = null;

	static public String detailNameChange = "DetailNameChange";

    public StsGroupBox displayBox;
	public StsGroupBox defaultBox;
	public StsGroupBox propertyBox;
//	JTextField detailNameField = new JTextField();
//	JLabel detailNameLabel = new JLabel();
    boolean editable = false;

    Class c = null;
	Vector actionListeners = null;
	Vector focusListeners = null;
    private transient Vector changeListeners;

    private StsObjectPanel(StsTreeObjectI treeObject, boolean editable) throws Exception
    {
        setObject(treeObject);
        this.editable = editable;

        displayFields = treeObject.getDisplayFields();
		defaultFields = treeObject.getDefaultFields();
        propertyFields = treeObject.getPropertyFields();
        setPanelInFields(displayFields);
        setPanelInFields(propertyFields);
		setPanelInFields(defaultFields);

        c = treeObject.getClass();

	    if(treeObject instanceof StsProject || treeObject instanceof StsClass)
		{
			displayBox = StsGroupBox.noInsets("Display Preferences");
            propertyBox = StsGroupBox.noInsets("Properties Preferences");
			defaultBox = StsGroupBox.noInsets("Default Preferences");
		}
		else
		{
			displayBox = StsGroupBox.noInsets("Display");
            propertyBox = StsGroupBox.noInsets("Properties");
			defaultBox = StsGroupBox.noInsets("Defaults");
		}

        jbInit();
    }

    static public StsObjectPanel constructor(StsTreeObjectI treeObject, boolean editable)
    {
        try
        {
            return new StsObjectPanel(treeObject, editable);
        }
        catch(Exception e)
        {
            StsException.outputException("StsObjectPanel.constructor() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

	static public StsObjectPanel checkConstruct(StsObjectPanel currentPanel, StsTreeObjectI treeObject, boolean editable)
	{
		try
		{
			if(currentPanel != null && currentPanel.hasPanelObject(treeObject))
				return currentPanel;
			else
				return new StsObjectPanel(treeObject, editable);
		}
		catch(Exception e)
		{
			StsException.outputException("StsObjectPanel.constructor() failed.",
				e, StsException.WARNING);
			return null;
		}
	}

    private StsObjectPanel()
    {
    }

    static public StsObjectPanel constructEmptyPanel()
    {
        try
        {
            return new StsObjectPanel();
        }
        catch(Exception e)
        {
            StsException.systemError("StsObjectPanel.constructEmptyPanel() failed.");
            return null;
        }
    }

    private void setPanelInFields(StsFieldBean[] fields)
    {
        if(fields == null) return;

        for(int n = 0; n < fields.length; n++)
        {
            if(fields[n] == null) continue;
            fields[n].setBeanPanel(this);
        }
    }

	private void jbInit() throws Exception
	{
		// pack the fields tight vertically
		gbc.insets = new Insets(0, 0, 0, 0);
		gbc.anchor = GridBagConstraints.WEST;
//		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.weightx = 0.0;
        gbc.fill = gbc.HORIZONTAL;

        if(displayFields != null)
        {
            addFields(displayBox, displayFields);
            add(displayBox);
            add(Box.createVerticalStrut(10));
        }

		if(defaultFields != null)
		{
			addFields(defaultBox, defaultFields);
			add(defaultBox);
			add(Box.createVerticalStrut(10));
		}
       if(propertyFields != null)
        {
            addFields(propertyBox, propertyFields);
			add(propertyBox);
            add(Box.createVerticalStrut(15));
        }
	}

    public void setObject(StsTreeObjectI object) { panelObjects = new StsTreeObjectI[] { object }; }
    public StsTreeObjectI getObject()
    {
        if(panelObjects == null) return null;
        else return panelObjects[0];
    }

    public void setObjects(StsTreeObjectI[] objects) { panelObjects = objects; }

    public Object getPanelObject()
    {
        if(panelObjects == null) return null;
        return panelObjects[0];
    }

	public boolean hasPanelObject(Object object)
	{
		if(panelObjects == null) return false;
		for(int n = 0; n < panelObjects.length; n++)
			if(panelObjects[n] == object) return true;
	    return false;
	}

    public Object[] getPanelObjects() { return panelObjects; }

    public Class getDisplayClass() { return null; }

	public synchronized void removeActionListener(ActionListener l)
	{
		if (actionListeners != null && actionListeners.contains(l))
		{
			Vector v = (Vector) actionListeners.clone();
			v.removeElement(l);
			actionListeners = v;
		}
	}

	public synchronized void addActionListener(ActionListener l)
	{
		Vector v = actionListeners == null ? new Vector(2) : (Vector) actionListeners.clone();
		if (!v.contains(l))
		{
			v.addElement(l);
			actionListeners = v;
		}
	}

	protected void fireActionPerformed(ActionEvent e)
	{
		if (actionListeners != null)
		{
			Vector listeners = actionListeners;
			int count = listeners.size();
			for (int i = 0; i < count; i++)
			  ((ActionListener) listeners.elementAt(i)).actionPerformed(e);
		}
	}

	public void actionPerformed(ActionEvent e)
	{
    	fireActionPerformed(e);
	}

	public void stateChanged(ChangeEvent e)
	{
//    	fireStateChanged(e);
	}

/*
    public StsFieldBean getDisplayField(Component c) { return getFieldForComponent(displayFields, c); }
    public StsFieldBean getPropertyField(Component c) { return getFieldForComponent(propertyFields, c); }

    public StsFieldBean getFieldForComponent(StsFieldBean[] fields, Component c)
    {
		if(fields == null) return null;

        for( int i=0; i<fields.length; i++ )
        {
            Component[] components = fields[i].getBeanComponents();
            if( components != null )
            {
                if( (components.length == 1 && components[0] == c) ||
                    (components.length == 2 && components[1] == c) )
                    return fields[i];
            }
        }

        return null;
    }
*/

    public StsFieldBean getDisplayField(String fieldName) { return getField(displayFields, fieldName); }
    public StsFieldBean getPropertyField(String fieldName) { return getField(propertyFields, fieldName); }

    public StsFieldBean getField(StsFieldBean[] fields, String otherFieldName)
    {

		if(fields == null) return null;

        for( int i=0; i<fields.length; i++ )
        {
            String fieldName = fields[i].getFieldName();
            if( fieldName != null && fieldName.equals(otherFieldName) )
                return fields[i];
        }
        return null;
    }

    public StsFieldBean[] getDisplayFields() { return displayFields; }
    public StsFieldBean[] getPropertyFields() { return null; }

    void addFields(StsGroupBox groupBox, StsFieldBean[] fields)
    {
		groupBox.gbc.anchor = GridBagConstraints.WEST;
//		groupBox.gbc.fill = GridBagConstraints.HORIZONTAL;
        if(fields == null) return;
		for( int i=0; i<fields.length; i++ )
		{
			groupBox.addBeanPanel(fields[i]);
//			Component component = fields[i].getObjectPanelComponent();
//			groupBox.add(component);
		}
		groupBox.gbc.anchor = GridBagConstraints.CENTER;
    }
/*
        GridBagLayout gb = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();

        groupBox.setLayout(gb);

        gbc.anchor = GridBagConstraints.WEST;
        for( int i=0; i<fields.length; i++ )
        {
            Component[] components = fields[i].getBeanComponents();
            components = removeNullComponents(components);

            gbc.gridy = i;
            int nComponents = components.length;
            if(nComponents == 0) continue;
            if(nComponents == 1)
                gbc.gridwidth = 2;
            else
                gbc.gridwidth = 1;

            groupBox.add(components[0], gbc);

            if(nComponents == 1) continue;

            groupBox.add(Box.createHorizontalStrut(10));

            gbc.gridwidth = GridBagConstraints.REMAINDER;
            gbc.fill = GridBagConstraints.HORIZONTAL;
            groupBox.add(components[1], gbc);
        }
    }

    private Component[] removeNullComponents(Component[] components)
    {
        if(components == null) return new Component[0];
        int nNotNull = 0;
        for(int n = 0; n < components.length; n++)
            if(components[n] != null) nNotNull++;

        if(nNotNull == components.length) return components;

        Component[] newComponents = new Component[nNotNull];
        int nn = 0;
        for(int n = 0; n < components.length; n++)
            if(components[n] != null) newComponents[nn++] = components[n];
        return newComponents;
    }
*/
    public void refreshProperties()
    {
        if( panelObjects == null ) return;
        Object panelObject = panelObjects[0];
        if(panelObject == null) return;

 //       refreshDetailName();
        displayFields = ((StsTreeObjectI)panelObject).getDisplayFields();
        if(displayFields != null)
        {
            for( int i=0; i<displayFields.length; i++ )
                displayFields[i].setValueFromPanelObject(panelObject);
            // displayBox.invalidate();
		}
        propertyFields = ((StsTreeObjectI)panelObject).getPropertyFields();
        if(propertyFields != null)
        {
            for( int i=0; i<propertyFields.length; i++ )
            {
                if(propertyFields[i] == null) continue;
                propertyFields[i].setValueFromPanelObject(panelObject);
            }
            // propertyBox.invalidate();
        }
        validate();
        // paintImmediately(0, 0, getWidth(), getHeight());
        this.repaint();
    }
/*
    public void refreshDetailName()
    {

        if( selected == null || selected.length == 0 )
        {
        	setDetailName("<none>");
//            enableAll(false);
        }
        else
        {
//        	enableAll(true);
            if(selected.length == 1 )
            {
	            refreshDetailName(getListItem(selected[0]));
            }
            else
            {
                setDetailName("<multiple>");
            }
        }
//        panel.setNameEditable( selected == null ? false : selected.length == 1);
	}

	public void refreshDetailName(Object obj)
    {
        try
        {
            String name = (String)getName.invoke(obj, null);
            setDetailName(name);
        }
        catch(Exception ex) { }
    }
*/
    Object getListItem(int index) { return getListItem(c, index); }
    Object getPropertyItem(int index) { return getListItem(index); }

    Object getListItem(Class listClass, int index)
    {
        if(panelObjects == null) return null;
        StsTreeObjectI panelObject = panelObjects[0];
        if(panelObject == null) return null;
        if(!(panelObject instanceof StsClass)) return null;
        StsClass stsClass = (StsClass)panelObject;
        return stsClass.getElement(index);
    }

    int getListSize() { return getListSize(c); }

    public int getListSize(Class listClass)
    {
        if(panelObjects == null) return 0;
        StsTreeObjectI panelObject = panelObjects[0];
        if(!(panelObject instanceof StsClass)) return 0;
        StsClass stsClass = (StsClass)panelObject;
        return stsClass.getSize();
	}

    public boolean isSelected(int index)
    {
    	if(selected == null) return false;
        for(int i=0; i<selected.length; i++ )
            if(selected[i] == index) return true;
        return false;
    }

    // doesn't work properly
    public void removeDisplayField(StsFieldBean fieldBean)
    {
        try
        {
            displayFields = removeField(displayFields, fieldBean);
            jbInit();
        }
        catch(Exception e)
        {
            StsException.systemError("StsObjectPanel.removeDisplayField() failed for fieldBean " + fieldBean.getName());
        }
    }

    // doesn't work properly
    private StsFieldBean[] removeField(StsFieldBean[] fieldBeans, StsFieldBean fieldBean)
    {
        return (StsFieldBean[])StsMath.arrayDeleteElement(fieldBeans, fieldBean);
    }

    // doesn't work properly
    public void addDisplayField(StsFieldBean fieldBean)
    {
        try
        {
            displayFields = addField(displayFields, fieldBean);
            jbInit();
        }
        catch(Exception e)
        {
            StsException.systemError("StsObjectPanel.removeDisplayField() failed for fieldBean " + fieldBean.getName());
        }
    }

    // doesn't work properly
    private StsFieldBean[] addField(StsFieldBean[] fieldBeans, StsFieldBean fieldBean)
    {
        return (StsFieldBean[])StsMath.arrayAddElement(fieldBeans, fieldBean);
    }


    public void setSelected(int[] indices)
    {
//    	list.setSelectedIndices(indices);
//        refreshDetailName(indices);
		selected = indices;
    }

/*
    public void updateFieldBean(String fieldName, Object value)
    {
        StsFieldBean field = getFieldBean(fieldName);
        field.setValue(value);
    }

    public void updateFieldBeans(Object object)
    {
        StsFieldBean[] fields = getPropertyFields();
        if(fields == null) return;
        for(int n = 0; n < fields.length; n++)
        {
            fields[n].setValue(object);
        }
    }

    public void updatePropertyField(String fieldName, Object value)
    {
        StsFieldBean field = getPropertyField(fieldName);
        field.setValue(value);
    }

    public void updatePropertyFields(Object object)
    {
        StsFieldBean[] fields = getPropertyFields();
        if(fields == null) return;
        for(int n = 0; n < fields.length; n++)
        {
            fields[n].setValue(object);
        }
    }

	public synchronized void removeFocusListener(FocusListener l)
	{
		if (focusListeners != null && focusListeners.contains(l))
		{
			Vector v = (Vector) focusListeners.clone();
			v.removeElement(l);
			focusListeners = v;
		}
	}

	public synchronized void addFocusListener(FocusListener l)
	{
		Vector v = focusListeners == null ? new Vector(2) : (Vector) focusListeners.clone();
		if (!v.contains(l))
		{
			v.addElement(l);
			focusListeners = v;
		}
	}

	protected void fireFocusGained(FocusEvent e)
	{
		if (focusListeners != null)
		{
			Vector listeners = focusListeners;
			int count = listeners.size();
			for (int i = 0; i < count; i++)
			  ((FocusListener) listeners.elementAt(i)).focusGained(e);
		}
	}

	protected void fireFocusLost(FocusEvent e)
	{
		if (focusListeners != null)
		{
			Vector listeners = focusListeners;
			int count = listeners.size();
			for (int i = 0; i < count; i++)
			  ((FocusListener) listeners.elementAt(i)).focusLost(e);
		}
	}

	public void focusGained(FocusEvent e) { fireFocusGained(e); }
	public void focusLost(FocusEvent e) { fireFocusLost(e); }

    public synchronized void removeChangeListener(ChangeListener l)
    {
        if(changeListeners != null && changeListeners.contains(l))
        {
            Vector v = (Vector) changeListeners.clone();
            v.removeElement(l);
            changeListeners = v;
        }
    }

    public synchronized void addChangeListener(ChangeListener l)
    {
        Vector v = changeListeners == null ? new Vector(2) : (Vector) changeListeners.clone();
        if(!v.contains(l))
        {
            v.addElement(l);
            changeListeners = v;
        }
    }

    protected void fireStateChanged(ChangeEvent e)
    {
        if(changeListeners != null)
        {
            Vector listeners = changeListeners;
            int count = listeners.size();
            for (int i = 0; i < count; i++)
            {
                ((ChangeListener) listeners.elementAt(i)).stateChanged(e);
            }
        }
    }
*/
    static public void main(String[] args)
    {
        try
        {
            ObjectPanelTest objectPanelTest = new ObjectPanelTest();
            StsObjectPanel panel = StsObjectPanel.constructor(objectPanelTest, true);
            StsToolkit.createDialog(panel, false);
//            objectPanelTest.printState();
//			panel = new StsObjectPanel(objectPanelTest, true);
//			StsToolkit.createDialog(panel);
//			objectPanelTest.printState();
        }
        catch(Exception e) { e.printStackTrace(); }
    }
}



