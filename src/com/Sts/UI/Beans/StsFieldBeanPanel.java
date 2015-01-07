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

public class StsFieldBeanPanel extends StsJPanel implements StsBeanPanelI
{
    protected Object[] panelObjects = null; // Object from which this panel get/set(s) values
//    protected StsFieldBean[] fieldBeans = null;
    private Component vertStrut5 = Box.createVerticalStrut(5);
    public StsFieldBeanPanel()
    {
//        initializeLayout();
//        setLayout();
    }

    /** default layout required for UI designers like JBuilder */
    private void setLayout()
    {
//        setLayout(xYLayout1);
        setSize(400, 300);
    }
/*
	public void initializeLayout()
	{
        try
        {
            setLayout(gb);
            gbc.insets = new Insets(2, 2, 2, 2);
            gbc.anchor = GridBagConstraints.WEST;
            gbc.gridwidth = 1;
            gbc.gridx = 0;
            gbc.gridy = 0;
            gbc.fill = GridBagConstraints.HORIZONTAL;
        }
        catch(Exception e)
        {
            StsException.systemError("StsFieldBeanPanel.initializeLayout() failed.");
        }
    }
*/

    /** when beans are added to panel, they are first given this beanPanel */
/*
    public void add(StsFieldBean fieldBean)
    {
        fieldBean.setBeanPanel(this);
        addField(fieldBean);
    }
*/
    public Component add(Component component)
    {
		checkAddFieldBeans(component);
        gbc.gridx = 0;
        add(component, gbc);
        return component;
    }

	public Component add(StsFieldBean bean)
	{
		bean.addToPanel(this);
		return bean;
	}

	private void checkAddFieldBeans(Component component)
	{
		if(component instanceof StsJPanel)
		{
			StsFieldBean[] newBeans = ((StsJPanel)component).getFieldBeans();
			if(newBeans == null) return;
			fieldBeans = (StsFieldBean[])StsMath.arrayAddArray(fieldBeans, newBeans, StsFieldBean.class);
		}
	}

    public void add(Component component, GridBagConstraints gbc)
    {
        super.add(component, gbc);
        gbc.gridy += 1;
        super.add(vertStrut5, gbc);
        gbc.gridy += 1;
    }
/*
    public void add(StsFieldBean field, GridBagConstraints gbc)
    {
        if(field == null) return;

        Component[] components = field.getNonNullComponents();
//        components = removeNullComponents(components);

        int nComponents = components.length;
        if(nComponents > 0)
        {
            super.add(components[0], gbc);
            gbc.gridy += 1;
        }
        if(nComponents > 1)
        {
            gbc.gridwidth = GridBagConstraints.REMAINDER;
            gbc.gridx += 1;
            gbc.anchor = GridBagConstraints.EAST;
            super.add(components[1], gbc);
            gbc.anchor = GridBagConstraints.WEST;
            gbc.gridx -= 1;
            gbc.gridwidth = 1;
       }
       gbc.gridy += 1;
	   super.add(vertStrut5, gbc);
	   gbc.gridy += 1;
	   fieldBeans = (StsFieldBean[])StsMath.arrayAddElement(fieldBeans, field, StsFieldBean.class);
    }

    public void remove(StsFieldBean field)
    {
        if(fieldBeans == null || fieldBeans.length == 0) return;
        fieldBeans = (StsFieldBean[])StsMath.arrayDeleteElement(fieldBeans, field);
        Component[] components = field.getBeanComponents();
        if(components == null) return;
        for(int n = 0; n < components.length; n++)
            if(components[n] != null) super.remove(components[n]);
    }

    private void addField(StsFieldBean field)
    {
        if(field == null) return;

        Component[] components = field.getNonNullComponents();
//        components = removeNullComponents(components);

        int nComponents = components.length;
        if(nComponents > 0)
        {
            gbc.gridwidth = 1;
            gbc.gridx = 0;
            super.add(components[0], gbc);
        }
        if(nComponents > 1)
        {
            gbc.gridwidth = GridBagConstraints.REMAINDER;
            gbc.gridx = 1;
            super.add(components[1], gbc);
            gbc.gridwidth = 1;
        }
        gbc.gridy += 1;
        super.add(vertStrut5, gbc);
          gbc.gridy += 1;

        fieldBeans = (StsFieldBean[])StsMath.arrayAddElement(fieldBeans, field, StsFieldBean.class);
    }
*/
    public StsFieldBean getBeanNamed(String name)
    {
        if(fieldBeans == null) return null;
        for(int n = 0; n < fieldBeans.length; n++)
            if(fieldBeans[n].getName().equals(name)) return fieldBeans[n];
        return null;
    }

    public void setPanelObject(Object object)
    {
        if(fieldBeans == null) return;
        panelObjects = new Object[] { object };
        for(int n = 0; n < fieldBeans.length; n++)
            fieldBeans[n].setBeanObject(object);
    }

    public Object getPanelObject()
    {
        if(panelObjects == null) return null;
        return panelObjects[0];
    }

    public void setPanelObjects(Object[] objects) { panelObjects = objects; }
    public Object[] getPanelObjects() { return panelObjects; }

    public void setEditable(boolean enabled)
    {
        if(fieldBeans == null) return;
        for(int n = 0; n < fieldBeans.length; n++)
            fieldBeans[n].setEditable(enabled);
    }
}
