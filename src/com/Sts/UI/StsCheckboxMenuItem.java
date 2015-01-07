
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.event.*;

/** Adds a convenient methods to set menuItems. */

public class StsCheckboxMenuItem extends JCheckBoxMenuItem implements Cloneable, ItemListener
{
	public StsCheckboxMenuItem()
	{
	}

    /** Sets the label and attaches an itemListener
      * @param label menuItem label
      * @param target instance object listening to to this menuItem
      * @param methodName String name of the instance method listening to this menuItem
      */
    public void setMenuItemListener(String label, Object target, String methodName)
    {
    	try
        {
    		setText(label);
            setName(label);
            addItemListener(new StsItemListener(target, methodName));
        }
        catch(Exception e)
        {
        	StsException.outputException(e, StsException.FATAL);
        }
    }

	/** For this menuItemListener, we pass the label as an additional argument.  This is used when the method being called
	 *  is a general toggle method and is not specific to the property being toggled.
	 */
	public void setToggleMenuItemListener(String label, Object target, String methodName)
	{
		try
		{
			setText(label);
			addItemListener( new StsItemListener(target, methodName, new Object[] {label}) );
		}
		catch(Exception e)
		{
			StsException.outputException("StsCheckboxMenuItem.setToggleMenuItemListener() failed.", e, StsException.FATAL);
		}
	}

    /** Sets the label and attaches an actionListener
      * @param label menuItem label
      * @param c class of the class/method listening to to this menuItem
      * @param methodName String name of the class/method listening to this menuItem
      */
    public void setMenuItemListener(String label, Class c, String methodName)
    {
    	try
        {
    		setText(label);
        	addItemListener(new StsItemListener(c, methodName));
        }
        catch(Exception e)
        {
        	StsException.outputException(e, StsException.FATAL);
        }
    }

    public StsCheckboxMenuItem copy()
    {
        try
        {
            StsCheckboxMenuItem menuItem = new StsCheckboxMenuItem();
            menuItem.setText(getText());
            menuItem.setSelected(this.isSelected());
            ItemListener[] itemListeners = getItemListeners();
            if(itemListeners != null)
            for(int n = 0; n < itemListeners.length; n++)
                menuItem.addItemListener(itemListeners[n]);
            menuItem.addItemListener(this);
            return menuItem;
        }
        catch (Exception e)
        {
            System.out.println("Exception in StsObject()\n" + e);
            return null;
        }
    }

    public void itemStateChanged(ItemEvent e)
    {
        boolean visible = (e.getStateChange() == ItemEvent.SELECTED);
        this.getModel().setSelected(visible);
    }
}
