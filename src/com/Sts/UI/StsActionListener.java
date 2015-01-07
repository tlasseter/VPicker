package com.Sts.UI;

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

import com.Sts.Reflect.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.awt.event.*;

public class StsActionListener implements ActionListener
{
  	static protected int clicks = 0; // Used for test routine only

  	protected StsMethod method;

  /** constructor
    * @param method StsMethod containing info for method invocation
    */
    public StsActionListener(StsMethod method)
    {
    	this.method = method;
    }

    public void actionPerformed(ActionEvent event)
    {
        try
        {
            method.invokeInstanceMethod();
        }
        catch (Exception e)
        {
			if(method != null)
				StsException.outputException( "method=" + method.getMethod().getName(), e, StsException.WARNING);
			else
				StsException.outputException( "method=null", e, StsException.WARNING);
        }
    }

    // A simple test program for the StsActionListener
    public static void main(String[] args) throws NoSuchMethodException
    {
        Frame f = new Frame("StsActionListener Test");      // Create window.
        f.setLayout(new FlowLayout());                      // Set layout manager.
        Button b1 = new Button("tick");                     // Create buttons.
        Button b2 = new Button("tock");
        Button b3 = new Button("Close Window");
        Button b4 = new Button("static");
        Button b5 = new Button("class");
        f.add(b1); f.add(b2); f.add(b3); f.add(b4); f.add(b5);  // Add them to window.

        // Specify what the buttons do.  Invoke a named method with
        // the StsActionListener object.

        b1.addActionListener(new StsActionListener(new StsMethod(b1, "setLabel", "tock")) );
        b1.addActionListener(new StsActionListener(new StsMethod(b2, "setLabel", "tick")) );
        b1.addActionListener(new StsActionListener(new StsMethod(b3, "hide", null)) );
        b2.addActionListener(new StsActionListener(new StsMethod(b1, "setLabel", "tick")) );
        b2.addActionListener(new StsActionListener(new StsMethod(b2, "setLabel", "tock")) );
        b2.addActionListener(new StsActionListener(new StsMethod(b3, "show", null)) );
        b3.addActionListener(new StsActionListener(new StsMethod(f, "dispose", null)) );

       	b5.addActionListener(new StsActionListener(new StsMethod("com.Sts.Types.StsActionListener", "staticMethod", new Object[]{b4,"dummy"} )) );

       	try
        {
        	Class c = Class.forName("com.Sts.Types.StsActionListener");
        	b4.addActionListener(new StsActionListener(new StsMethod(c, "staticMethod", b5)) );
        }
        catch (Exception e) { System.err.println("StsActionListener: " + e); }

        f.pack();                                             // Set window size.
        f.setVisible(true);                                             // And pop it up.
    }

  	static void staticMethod(Button b)
  	{
  		clicks++;
    	b.setLabel(Integer.toString(clicks));
  	}

    static void staticMethod(Button b, String name) // Tests multiple arguments
  	{
  		clicks++;
    	b.setLabel(Integer.toString(clicks));
  	}
}

