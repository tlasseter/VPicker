//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC


//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.MVC.*;
import com.Sts.Reflect.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class StsToggleButton extends JToggleButton
{
    StsMethod toggleMethod = null;
    Icon selectedIcon = null, deselectedIcon = null;
    String selectedName = null, deselectedName = null;

    /** So as to not cover the button, set the tooltip ABOVE or BELOW.  For normal
     *  tooltip display, set topPosition to CENTER
     */
    private byte tipPosition = BELOW;

    public static final byte ABOVE = 0;
    public static final byte BELOW = 1;
    public static final byte CENTER = 2;
    private static Font font = UIManager.getFont("ToolTip.font");

    // Build a custom help cursor
    static
    {
        try
        {
            Class types[] = {Image.class, Point.class, String.class};
            Toolkit.class.getMethod("createCustomCursor", types);
        }
        catch (Exception e)
        {System.out.println("Unable to create custom help cursor.");
        }
    }

	static final int defaultWidth = 20;
	static final int defaultHeight = 20;

    static final boolean debug = false;

    public StsToggleButton()
    {
    }

    /** constructor for instance listener instance.methodName with args list */
    public StsToggleButton(String name, String tip)
    {
        super();

        try
        {
            labelButton(name);
            setName(name);
            setToolTipText(tip);
            setMargin(new Insets(0, 0, 0, 0));
        }
        catch (Exception e)
        {
            StsException.systemError("StsToggleButton.constructor(name, tip) failed.");
        }
    }

    /** constructor for class listener actionClass.newAction with list of args */
    public StsToggleButton(String name, String tip, Class actionClass, Object[] args)
    {
        this(name, tip, actionClass, null, args);
    }

    /** constructor for class listener actionClass.method with single arg */
    public StsToggleButton(String name, String tip, Class actionClass, String methodName, Object arg)
    {
        this(name, tip, actionClass, methodName, (arg == null) ? null : new Object[]
             {arg});
    }

    /** constructor for class listener actionClass.method with args list */
    public StsToggleButton(String name, String tip, Class actionClass, String methodName, Object[] args)
    {
        super();

        try
        {
            labelButton(name);
            setName(name);
            setToolTipText(tip);
            setMargin(new Insets(0, 0, 0, 0));
            addActionListener(actionClass, methodName, args);
        }
        catch (Exception e)
        {
            StsException.systemError("StsToggleButton.constructor(class.method(args)) failed.");
        }
    }

    /** constructor for instance listener instance.methodName with no args */
    public StsToggleButton(String name, String tip, Object instance, String methodName)
    {
        this(name, tip, instance, methodName, (Object[])null);
    }

    /** constructor for instance listener instance.methodName with single arg */
    public StsToggleButton(String name, String tip, Object instance, String methodName, Object arg)
    {
        this(name, tip, instance, methodName, (arg == null) ? null : new Object[]{arg} );
    }


    /** constructor for instance listener instance.methodName with single arg */
    public StsToggleButton(String name, String tip, Object instance, String selectMethodName, String deselectMethodName)
    {
        this(name, tip, instance, selectMethodName, deselectMethodName, null);
    }

    /** constructor for instance listener instance.methodName with single arg */
    public StsToggleButton(String name, String tip, Object instance, String selectMethodName, String deselectMethodName, Object arg)
    {
        this(name, tip, instance, selectMethodName, deselectMethodName, (arg == null) ? null : new Object[] {arg});
    }

    /** constructor for instance listener instance.methodName with args list */
    public StsToggleButton(String name, String tip, Object instance, String methodName, Object[] args)
    {
        super();

        try
        {
            labelButton(name);
            setName(name);
            setToolTipText(tip);
            setMargin(new Insets(0, 0, 0, 0));
            StsMethod m = new StsMethod(instance, methodName, args);
            addActionListener(new StsActionListener(m));
        }
        catch (Exception e)
        {
            StsException.systemError("StsToggleButton.constructor(instance.method(args)) failed.");
        }
    }

    /** constructor for instance listener instance.methodName with args list */
    public StsToggleButton(String name, String tip, Object instance, String selectMethodName, String deselectMethodName, Object[] args)
    {
        super();

        try
        {
            labelButton(name);
            setName(name);
            setToolTipText(tip);
            setMargin(new Insets(0, 0, 0, 0));
            StsMethod selectMethod = new StsMethod(instance, selectMethodName, args);
            StsMethod deselectMethod = new StsMethod(instance, deselectMethodName, args);
            addActionListener(new StsToggleActionListener(this, selectMethod, deselectMethod));
        }
        catch (Exception e)
        {
            StsException.systemError("StsToggleButton.constructor(instance.method(args)) failed.");
        }
    }

    /** constructor for actionManger.newAction with no args */
    public StsToggleButton(String name, String tip, StsActionManager actionManager, Class actionClass)
    {
        super();
        initialize(name, tip, actionManager, actionClass, null);
    }

    /** constructor for actionManger.newAction with single arg */
    public StsToggleButton(String name, String tip, StsActionManager actionManager, Class actionClass, Object arg)
    {
        super();
        initialize(name, tip, actionManager, actionClass, (arg == null) ? null : new Object[]
                   {arg});
    }

    /** constructor for actionManager.newAction with args list */
    public StsToggleButton(String name, String tip, StsActionManager actionManager, Class actionClass,
                           Object[] args)
    {
        super();
        initialize(name, tip, actionManager, actionClass, args);
    }

    public void initialize(String name, String tip, StsActionManager actionManager, Class actionClass, Object arg)
    {
        initialize(name, tip, actionManager, actionClass, (arg == null) ? null : new Object[] {arg});
    }

    public void initialize(String name, String tip, StsActionManager actionManager, Class actionClass)
    {
        initialize(name, tip, actionManager, actionClass, null);
    }

    public void initialize(String name, String tip, StsActionManager actionManager, Class actionClass, Object[] args)
    {
        try
        {
            if (actionManager == null)
            {
                StsException.systemError("StsToggleButton.constructor(actionManger..) failed." +
                                         " ActionManager cannot be null.");
                return;
            }
            labelButton(name);
            setName(name);
            setToolTipText(tip);
            setMargin(new Insets(0, 0, 0, 0));

            Object[] selectArguments = new Object[] {actionClass, (args == null) ? new Object[0] : args};
            StsMethod selectMethod = new StsMethod(actionManager, "startAction", selectArguments);

            Object[] deselectArguments = new Object[] {actionClass};
            StsMethod deselectMethod = new StsMethod(actionManager, "endAction", deselectArguments);

            addActionListener(new StsToggleActionListener(this, selectMethod, deselectMethod));

            actionManager.addActionWorkflow(actionClass);
        }
        catch (Exception e)
        {
            StsException.systemError("StsToggleButton.constructor(class.method(args)) failed.");
        }
    }

    private void addActionListener(Class actionClass, String methodName, Object[] args)
    {

        StsMethod m = new StsMethod(actionClass, methodName, args);
        super.addActionListener(new StsActionListener(m));
    }

    public void addActionListener(Object instance, String methodName, Object[] args)
    {
        StsMethod m = new StsMethod(instance, methodName, args);
        addActionListener(new StsActionListener(m));
    }

    public void addActionListener(Object instance, String methodName, Object arg)
    {
        StsMethod m = new StsMethod(instance, methodName, new Object[] {arg});
        addActionListener(new StsActionListener(m));
    }

    private void labelButton(String name)
    {
        Icon icon = StsIcon.createIcon(name + ".gif");
        if (icon != null)
        {
            setIcon(icon);
        }
        else
        {
            setText(name);
        }
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

    public void addIcons(String selectedName, String deselectedName)
    {
        setText(null);
        selectedIcon = StsIcon.createIcon(selectedName + ".gif");
        if (selectedIcon == null)
        {
            this.selectedName = selectedName;
        }

        deselectedIcon = StsIcon.createIcon(deselectedName + ".gif");
        if (deselectedIcon != null)
        {
            setIcon(deselectedIcon);
        }
        else
        {
            this.deselectedName = deselectedName;
            setText(deselectedName);
        }

        setSize(20, 20);

        addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                JToggleButton button = (JToggleButton) e.getSource();
                setIconOrName(button.isSelected());
            }
        });
    }

    private void setIconOrName(boolean selected)
    {
        if (selected)
        {
            if (selectedIcon != null)
            {
                setIcon(selectedIcon);
            }
            else if (selectedName != null)
            {
                setText(selectedName);
            }
        }
        else
        {
            if (deselectedIcon != null)
            {
                setIcon(deselectedIcon);
            }
            else if (deselectedName != null)
            {
                setText(deselectedName);
            }
        }
    }

    /*
        public Point getToolTipLocation(MouseEvent e)
        {
            Point location = e.getPoint();
            location.y += 16;
            location.x += 14;
            if( font != null )
            {
                int toolTipHeight = font.getSize() + 6;
             int height = location.y + toolTipHeight;
             if(height > getHeight()) location.y = getHeight() - toolTipHeight;
            }

            if( location.x > getWidth() ) location.x = getWidth() - 2;
         return location;
        }
     */
    public Point getToolTipLocation(MouseEvent e)
    {
        Point location = e.getPoint();

        if (tipPosition == ABOVE && font != null)
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

    public void fireActionPerformed()
    {
        super.fireActionPerformed(new ActionEvent(this, 0, this.getName()));
    }

    public void fireActionPerformed(ActionEvent event)
    {
        super.fireActionPerformed(event);
    }

    public void setSelected(boolean b)
    {
        if (debug)
        {
            boolean isSelected = isSelected();
            System.out.println("setSelected(" + b + ") called. Currently isSelected: " + isSelected);
        }
        setIconOrName(b);
        super.setSelected(b);
		// ?? fireActionPerformed(); // jbw
    }

    static StsToggleButton toggleButton;
    static StsToggleButton toggleButton1, toggleButton2;

    static public void main(String[] args)
    {
        try
        {
            javax.swing.UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
            System.out.println("Using Windows look-and-feel...");
            JFrame frame = new JFrame("StsToggleButton test");
            frame.getContentPane().setLayout(new BorderLayout());
            JButton button = new JButton("Change Toggle");
            StsToggleButtonTest test = new StsToggleButtonTest();

            toggleButton1 = new StsToggleButton("Toggle1", "do toggle1", test, "toggleOn", "toggleOff");
            toggleButton1.addIcons("xSelected", "xDeselected");
            Dimension size = new Dimension(12, 20);
            toggleButton1.setPreferredSize(size);
            toggleButton1.setMaximumSize(size);
 //           toggleButton2 = new StsToggleButton("Toggle2", "do toggle2", test, "toggleOn", "toggleOff");
			toggleButton1.setBorderPainted(true);
            StsJPanel panel = StsJPanel.addInsets();
            panel.gbc.fill = GridBagConstraints.NONE;
            panel.add(toggleButton1);
 //           ButtonGroup buttonGroup = new ButtonGroup();
 //           buttonGroup.add(toggleButton1);
 //           buttonGroup.add(toggleButton2);
            StsToolkit.createDialog(panel);
            frame.getContentPane().add(toggleButton1, BorderLayout.WEST);
 //           frame.getContentPane().add(toggleButton2, BorderLayout.EAST);
            frame.setSize(100, 50);
            frame.setVisible(true);
/*
            toggleButton1.setSelected(true);
            try { Thread.sleep(2000); } catch(Exception e) { }
            toggleButton2.setSelected(true);
            try { Thread.sleep(2000); } catch(Exception e) { }
            toggleButton1.setSelected(true);
            try { Thread.sleep(2000); } catch(Exception e) { }
*/
        }
        catch (Exception e)
        {
        }
    }
/*

    static public void main(String[] args)
    {
        try
        {
            javax.swing.UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
            System.out.println("Using Windows look-and-feel...");
            JFrame frame = new JFrame("StsToggleButton test");
            frame.getContentPane().setLayout(new BorderLayout());
            JButton button = new JButton("Change Toggle");
            StsToggleButtonTest test = new StsToggleButtonTest();

            toggleButton = new StsToggleButton("seismic", "toggle visiblity of seismic", test, "toggleOn", "toggleOff", null);
//        toggleButton = new StsToggleButton("Toggle", "This will do something");

            button.addActionListener
                (
                new ActionListener()
            {
                public void actionPerformed(ActionEvent e)
                {
                    boolean isPressed = toggleButton.getModel().isPressed();
//                    toggleButton.getModel().setPressed(!isPressed);
                    ActionEvent toggleEvent = new ActionEvent(new Object(), 0, "selected");
                    toggleButton.fireActionPerformed(toggleEvent);
                }
            }
            );

            frame.getContentPane().add(button, BorderLayout.WEST);
            frame.getContentPane().add(toggleButton, BorderLayout.EAST);
            frame.setSize(100, 50);
            frame.setVisible(true);
        }
        catch (Exception e)
        {
        }
    }
*/
}

class StsToggleButtonTest
{
    public StsToggleButtonTest()
    {
    }

    public void methodA()
    {
        System.out.println("Test.methodA() called.");
    }

    public void toggleOn()
    {
        System.out.println("Toggle on");
    }

    public void toggleOff()
    {
        System.out.println("Toggle off");
    }
}
