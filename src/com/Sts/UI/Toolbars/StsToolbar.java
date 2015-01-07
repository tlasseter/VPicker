
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Toolbars;

import com.Sts.DB.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.io.*;

//import com.Sts.Actions.Wizards.Collaboration.*;
/**
    The StsToolbar class is a generic S2S toolbar and has methods
    to set button icons, tool tips, and actions they perform.
    For a given StsToolbarPanel on an StsWin3dBase, add ALL possible toolbars
    to the panel.  As toolbars are turned on/off, their visiblity is toggled,
    NOT their existence.  This way a set of checkboxMenuItems is maintained which
    parallels the set of toolbars; toggle checkbox toggles the toolbar;
    closing the toolbar unchecks it on the menu list.
*/
public class StsToolbar extends JToolBar implements StsSerializable
{
    transient private StsWin3dBase parentWindow = null;
    protected boolean isVisible; //separate from JComponent.visible

    public StsToolbar()
	{
        initialize();
    }

    /** constructor called for a new toolbar*/
    public StsToolbar(String name)
	{
        this.setName(name);
        isVisible = getDefaultVisible();
        initialize();
    }

    public boolean getDefaultFloatable() { return true; }
    public boolean getDefaultVisible() { return true; }

    public void initialize()
    {
        setFloatable(getDefaultFloatable());
        setBorder(BorderFactory.createEmptyBorder());
        setToolTipText(getName());
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        return true;
    }
    /** displayed is redundant with visible in superClass, but is used for persistence of the toolbar visiblity state */
    public void setVisible(boolean visible)
    {
        super.setVisible(visible);
        this.isVisible = visible;

        //parentWindow.model.setProperty("isVisible", isVisible);  // Needs to be in the database for the viewer so toolbar state is included in archive
    }

    public boolean isVisible() { return isVisible; }

    /** convenience method to set the minimum size from the preferred size */
    protected void setMinimumSize() { setMinimumSize(getPreferredSize()); }

    /** enable/disable a button by its name */
    public boolean setButtonEnabled(Component buttonComponent, boolean enabled)
    {
        Component[] components = this.getComponents();

        for (int i = 0; i < components.length; i++)
        {
            if(components[i] == buttonComponent)
            {
                AbstractButton button = (AbstractButton)buttonComponent;
                button.setSelected(enabled);
                return true;
            }
        }
        return false;
    }

    /** Set visibility by button name */
    public boolean setButtonVisibility(Component buttonComponent, boolean visible)
    {
        Component[] components = this.getComponents();

        for (int i = 0; i < components.length; i++)
        {
            if(components[i] == buttonComponent)
            {
                AbstractButton button = (AbstractButton)buttonComponent;
                button.setVisible(visible);
                return true;
            }
        }
        return false;
    }

    public AbstractButton getButtonNamed(String buttonName)
    {
        Component[] components = this.getComponents();
        for (int i = 0; i < components.length; i++)
        {
            if(components[i] instanceof AbstractButton)
            {
                AbstractButton button = (AbstractButton)components[i];
                if(button.getName().equals(buttonName)) return button;
            }
        }
        return null;
    }

    static public StsToolbar getToolbarNamed(StsToolbar[] toolbars, String name)
    {
        if(toolbars.length == 0) return null;
        for(int n = 0; n < toolbars.length; n++)
            if(toolbars[n].getName().equals(name))
                return toolbars[n];
        return null;
    }

    public Component getComponentNamed(String name)
    {
        Component[] components = this.getComponents();
		if (name == null) System.out.println("getComponentNamed null");
        for (int i = 0; i < components.length; i++)
		{
			if (components[i] != null)
				if (components[i].getName() != null)
			       if(components[i].getName().equals(name))return components[i];
		}
        return null;
    }

    public void buttonTest(Integer i)
    {
        System.out.println("Button " + i.toString() + " pressed.");
    }

	public void toggleTest()
	{
		System.out.println("Toggle Button pressed.");
    }

    public void addComponent(Component component)
    {
        add(component);
    }

    public void replaceComponent(Component oldComponent, Component newComponent)
    {
        int index = getComponentIndex(oldComponent);
        remove(oldComponent);
        add(newComponent, index);
    }

    public boolean removeComponent(Component component)
    {
        Component[] components = this.getComponents();
        for (int i = 0; i < components.length; i++)
        {
            if(components[i] == component) components[i].setEnabled(false);
            return true;
        }
        return false;
    }

	public boolean selectToolbarItem(String itemName, boolean selected)
	{
		Component[] components = this.getComponents();
		for (int i = 0; i < components.length; i++)
		{
			if (components[i].getName().equals(itemName))
			{
				if (components[i] instanceof AbstractButton)
				{
					AbstractButton button = (AbstractButton) components[i];
					button.setSelected(selected);
					return true;
				}
			}
		}
		return false;
	}

	public void addCloseIcon(StsWin3dBase window)
	{
        parentWindow = window;
        Icon closeIcon = StsIcon.createIcon("toolbarClose.gif");
		add(new StsButton("close", closeIcon, "Close this toolbar.", this, "close"));
	}

	public void close()
	{
		parentWindow.closeToolbar(this);
	}

	public void writeObject(StsDBOutputStream out) throws IllegalAccessException, IOException
    {
        out.writeUTF(getName());
        out.writeUTF(getClass().getName());
        out.writeBoolean(isVisible());
    }

	public void readObject(StsDBInputStream in) throws IllegalAccessException, IOException
	{
        String name = in.readUTF();
        String classname = in.readUTF();
        boolean isVisible = in.readBoolean();
    }

    public boolean validateState() { return true; }
    public void reconfigure() {}

    static public void main(String[] args)
    {
        try
        {
			javax.swing.UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
			JDialog dialog = new JDialog((Frame)null, "Button Panel Test", true);
            StsToolbar toolbar = new StsToolbar("Toolbar Test");
            dialog.getContentPane().add(toolbar);
            toolbar.add(new StsButton("Button1", "This is button 1", toolbar, "buttonTest", new Integer(1)));
			StsToggleButton toggleButton = new StsToggleButton("Toggle", "This is a toggle button", toolbar, "toggleTest");
			Border buttonBorder = javax.swing.plaf.metal.MetalBorders.getButtonBorder();
			toggleButton.setBorder(buttonBorder);
//			LookAndFeel.installBorder(toggleButton, "Button.border");
//			toggleButton.setBorder(BorderFactory.createRaisedBevelBorder());
//			toggleButton.setBorderPainted(true);
			toolbar.add(toggleButton);
//			toolbar.addCloseIcon(dialog);
            dialog.pack();
            dialog.setVisible(true);
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

    public boolean forViewOnly() { return false; }
}
