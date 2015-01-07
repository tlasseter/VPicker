

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.Help.*;
import com.Sts.MVC.*;
import com.Sts.Reflect.*;
import com.Sts.UI.Icons.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class StsButton extends JButton
{
    /** So as to not cover the button, set the tooltip ABOVE or BELOW.  For normal
     *  tooltip display, set topPosition to CENTER
     */
    private byte tipPosition = BELOW;

    public static final byte ABOVE = 0;
    public static final byte BELOW = 1;
    public static final byte CENTER = 2;

   	private static Font font = UIManager.getFont("ToolTip.font");
    private static final String HELP_TARGET = "ui.toolbars.common";

    // Build a custom help cursor
    static
    {
    	try
        {
    	    Class types[] = {Image.class, Point.class, String.class};
	        Toolkit.class.getMethod("createCustomCursor", types);
    	}
        catch (Exception e) { System.out.println("Unable to create custom help cursor."); }
    }

	static final int defaultWidth = 20;
	static final int defaultHeight = 20;

	public StsButton()
	{
	}

    /** constructor for class listener actionClass.startAction with no args */
	public StsButton(String name, String tip, Class actionClass)
	{
        this(name, tip, actionClass, null, new Object[0] );
    }

    /** constructor for class listener actionClass.startAction with list of args */
	public StsButton(String name, String tip, Class actionClass, Object[] args)
	{
        this(name, tip, actionClass, null, args);
    }

    /** constructor for class listener actionClass.method with single arg */
	public StsButton(String name, String tip, Class actionClass, String methodName, Object arg)
	{
        this(name, tip, actionClass, methodName, (arg==null) ? null : new Object[]{arg});
    }

    /** constructor for class listener actionClass.method with args list */
	public StsButton(String name, String tip, Class actionClass, String methodName, Object[] args)
	{
        super();
        initialize(name, tip, actionClass, methodName, args);
	}

    /** constructor for instance listener instance.methodName with no args */
	public StsButton(String name, String tip, Object instance, String methodName)
	{
        this(name, tip, instance, methodName, null);
    }

    /** constructor for instance listener instance.methodName with no args */
	public StsButton(String name, Object instance, String methodName)
	{
        this(name, null, instance, methodName, null);
    }

    /** constructor for instance listener instance.methodName with single arg */
	public StsButton(String name, String tip, Object instance, String methodName, Object arg)
	{
        this(name, tip, instance, methodName, (arg==null) ? null : new Object[]{arg});
    }

    /** constructor for instance listener instance.methodName with args list */
	public StsButton(String name, String tip, Object instance, String methodName, Object[] args)
	{
        super();
        initialize(name, tip, instance, methodName, args);
    }

    /** constructor for actionManger.startAction with no args */
	public StsButton(String name, String tip, StsActionManager actionManager, Class actionClass)
	{
        this(name, tip, actionManager, actionClass, null);
    }

    /** constructor for actionManger.startAction with single arg */
	public StsButton(String name, String tip, StsActionManager actionManager, Class actionClass, Object arg)
	{
        this(name, tip, actionManager, actionClass, (arg==null) ? null : new Object[]{arg});
    }

    /** constructor for actionManger.startAction with args list */
	public StsButton(String name, String tip, StsActionManager actionManager, Class actionClass, Object[] args)
	{
        super();

    	try
        {
            if(actionManager == null)
            {
                StsException.systemError("StsButton.constructor(actionManger..) failed. ActionManager cannot be null.");
                return;
            }
		    labelButton(name);
            setName(name);
		    setToolTipText(tip);
		    setMargin(new Insets(0, 0, 0, 0));

            StsMethod m = null;
            Object[] arguments = new Object[]{actionClass, (args==null) ? new Object[0] : args};
            addActionListener(actionManager, "startAction", arguments);
        }
        catch(Exception e)
        {
        	StsException.systemError("StsButton.constructor(class.method(args)) failed.");
        }
	}

	/** constructor for instance listener instance.methodName with args list and icon label */
	public StsButton(String name, Icon icon, String tip, Object instance, String methodName, Object[] args)
	{
		super();
		initialize(name, icon, tip, instance, methodName, args);
    }

	/** constructor for instance listener instance.methodName with single arg and icon label */
	public StsButton(String name, Icon icon, String tip, Object instance, String methodName, Object arg)
	{
		super();
		initialize(name, icon, tip, instance, methodName, new Object[] { arg });
    }

	/** constructor for instance listener instance.methodName with no arg and icon label */
	public StsButton(String name, Icon icon, String tip, Object instance, String methodName)
	{
		this(name, icon, tip, instance, methodName, null);
    }

	public void initialize(String name, Icon icon, String tip, Object instance, String methodName, Object[] args)
	{
		try
		{
			setIcon(icon);
			setName(name);
			setToolTipText(tip);
			setMargin(new Insets(0, 0, 0, 0));
			StsMethod m = new StsMethod(instance, methodName, args);
			addActionListener(new StsActionListener(m));
		}
		catch(Exception e)
		{
			StsException.systemError("StsButton.constructor(instance.method(args)) failed.");
		}
    }

	public StsButton(String name, Icon icon, String tip, StsMethod method)
	{
		try
		{
			setIcon(icon);
			setName(name);
			setToolTipText(tip);
			setMargin(new Insets(0, 0, 0, 0));
			addActionListener(new StsActionListener(method));
		}
		catch(Exception e)
		{
			StsException.systemError("StsButton.constructor(instance.method(args)) failed.");
		}
    }

    public void initialize(String name, String tip, Class actionClass, String methodName, Object[] args)
    {

    	try
        {
            labelButton(name);
            setName(name);
            setToolTipText(tip);
            setMargin(new Insets(0, 0, 0, 0));
            addActionListener(actionClass, methodName, args);
        }
        catch(Exception e)
        {
        	StsException.systemError("StsButton.constructor(class.method(args)) failed.");
        }
    }

    public void initialize(String name, Object instance, String methodName)
    {
        initialize(name, null, instance, null);
    }

    public void initialize(String name, String tip, Object instance, String methodName)
    {
        initialize(name, tip, instance, null);
    }
    
    public void initialize(String name, String tip, Object instance, String methodName, Object[] args)
    {
    	try
        {
		    labelButton(name);
            setName(name);
		    setToolTipText(tip);
		    setMargin(new Insets(0, 0, 0, 0));
			StsMethod m = new StsMethod(instance, methodName, args);
        	addActionListener(new StsActionListener(m));
        }
        catch(Exception e)
        {
        	StsException.systemError("StsButton.constructor(instance.method(args)) failed.");
        }
    }

	public void addActionListener(Object instance, String methodName)
	{
        StsMethod m = new StsMethod(instance, methodName, null);
        addActionListener(new StsActionListener(m));
    }

	public void addActionListener(Object instance, String methodName, Object arg)
	{
        StsMethod m = new StsMethod(instance, methodName, new Object[] { arg } );
        addActionListener(new StsActionListener(m));
    }

	public void addActionListener(Object instance, String methodName, Object[] args)
	{
        StsMethod m = new StsMethod(instance, methodName, args);
        addActionListener(new StsActionListener(m));
    }

    private void addActionListener(Class actionClass, String methodName, Object[] args)
    {

        StsMethod m = new StsMethod(actionClass, methodName, args);
        super.addActionListener(new StsActionListener(m));
    }

    private void labelButton(String name)
    {
		Icon icon = StsIcon.createIcon(name + ".gif");
		if (icon != null)setIcon(icon);
		else setText(name);
		setSize();
	}

	public void setSize()
	{
		super.setSize(defaultWidth, defaultHeight);
	}

	public void setSize(int width, int height)
	{
		super.setSize(width, height);
	}


    static public StsButton constructHelpButton(Component parent)
    {
    	try
        {
		    StsButton helpButton = new StsButton();
            helpButton.labelButton("help");
		    helpButton.setToolTipText("Context-sensitive help");
            helpButton.setSize(20, 20);
		    helpButton.setMargin(new Insets(0, 0, 0, 0));
            HelpManager.setContextSensitiveHelp(helpButton, HelpManager.GENERAL, HELP_TARGET, parent);
	        HelpManager.setTrackingHelpActionListener(helpButton, HelpManager.GENERAL, parent);
            return helpButton;
        }
        catch(Exception e)
        {
        	StsException.systemError("StsButton.constructHelpButton() failed.");
            return null;
        }
    }

    public Point getToolTipLocation(MouseEvent e)
    {
        Point location = e.getPoint();

        if(tipPosition == ABOVE && font != null)
        {
            int toolTipHeight = font.getSize() + 6;
            int buttonHeight = getHeight();
            location.y -= buttonHeight + toolTipHeight;
        }
        else // tipPosition == BELOW
        {
            int buttonHeight = getHeight();
            location.y = buttonHeight;
        }
        return location;
    }

    public void setTipPositionAbove() { tipPosition = ABOVE; }
    public void setTipPositionBelow() { tipPosition = BELOW; }

    static public void main(String[] args)
    {
        JPanel panel = new JPanel();
        panel.setLayout(new BorderLayout());
//        panel.setSize(200, 200);
        StsButton button = new StsButton("button", "Button tip: push it.", null, "methodName", null);
        button.setPreferredSize(new Dimension(40, 15));
        button.setTipPositionAbove();
        panel.add(button, BorderLayout.CENTER);
        StsToolkit.createDialog(panel, true, 200, 200);
    }
}

