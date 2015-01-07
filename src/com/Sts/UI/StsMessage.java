
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.MVC.*;
import com.Sts.Reflect.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/** Use JOptionPane to display a message of type WARNING, INFO, or FLASH.  A
 *  FLASH message stays up for a few seconds and automatically disappears.  The
 *  other two types have an OK button to dismiss.  The message shouldn't have
 *  more than 32 characters or so per line or it will be truncated by JOptionPane:
 *  Insert a "\n" to make a new line.
 */
public class StsMessage implements Runnable
{
 	public static final int WARNING = 0;
    public static final int INFO = 1;
    public static final int FATAL = 2;
    public static final int ERROR = 3;

	private Component component;
	private int messageLevel;
	private String messageString;

    static final int noTypes = 5;

    static final String[] titles = new String[] { "Warning!!", "Information", "Fatal Error", "System Error" };



	public void run()
	{

		switch(messageLevel)
			{
				case WARNING:

					JOptionPane.showMessageDialog(component, messageString, titles[messageLevel],
												  JOptionPane.WARNING_MESSAGE);
					break;
				case INFO:

					JOptionPane.showMessageDialog(component, messageString, titles[messageLevel],
												  JOptionPane.INFORMATION_MESSAGE);
					break;
				case FATAL:

					String msg = messageString + "\nFATAL! Terminating application!";
					JOptionPane.showMessageDialog(component, msg, titles[messageLevel],
												  JOptionPane.ERROR_MESSAGE);
					System.exit(0);
					break;
				case ERROR:

					JOptionPane.showMessageDialog(component, messageString, titles[messageLevel],
												  JOptionPane.ERROR_MESSAGE);
					break;


			}


	}

	public StsMessage(Component component, int messageLevel, String messageString)
	{
        StsMessageFiles.errorMessage(toString(messageLevel) + "  " + messageString);
		this.component = component;
		this.messageLevel = messageLevel;
		this.messageString = messageString;
        StsToolkit.runLaterOnEventThread(this);
	}

	public StsMessage(Component component, int messageLevel, String messageString, boolean modal)
	{
		StsMessageFiles.errorMessage(toString(messageLevel) + "  " + messageString);
		this.component = component;
		this.messageLevel = messageLevel;
		this.messageString = messageString;
        if(!modal)
            StsToolkit.runLaterOnEventThread(this);
        else
			run();
	}

    /** Not currently used: alternate constructor with an additional message method
     *  useful if we wish StsMessage to route the message.  Currently it is the
     *  responsibility of the calling object to print the message anywhere else.
     */
	public StsMessage(Component component, String messageString, String titleString,
                      int messageIconType, StsMethod altMessageMethod)
	{

        try
        {
            JOptionPane.showMessageDialog(component, messageString, titleString, messageIconType);

            if(altMessageMethod != null)
                altMessageMethod.getMethod().invoke(altMessageMethod.getInstance(),
                                                    altMessageMethod.getMethodArgs());
        }
        catch(Exception e)
        {
            StsException.outputException("StsMessage.constructor() failed.",
                e, StsException.WARNING);
        }
	}

    static public String toString(int level)
    {
        switch (level)
        {
            case WARNING:   return titles[WARNING];
            case INFO:      return titles[INFO];
            case FATAL:     return titles[FATAL];
            case ERROR:     return titles[ERROR];
        }
        return "Unknown Message";
    }

    public static void main(String[] args)
    {
        try
        {
            UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
            //UIManager.setLookAndFeel(new com.sun.java.swing.plaf.motif.MotifLookAndFeel());
            //UIManager.setLookAndFeel(new javax.swing.plaf.metal.MetalLookAndFeel());
        }
        catch (Exception e)
        {
        }

        Frame frame = new Frame("Test Frame");
        frame.setSize(new Dimension(600, 400));
        StsToolkit.centerComponentOnScreen(frame);
        frame.setVisible(true);

        String message;

        message = new String("Info test: \nxxxxxxxxxxxx\nxxxxxxxx.");
        new StsMessage(frame, INFO, message);

        message = new String("Warning test: \nxxxxxxxxxxxx\nxxxxxxxx.");
        new StsMessage(frame, WARNING, message);

//        StsMethod altMessageMethod = new StsMethod(StsMessage.class, "printMessage", warningMessage);
    }

    static public void printMessage(String message)
    {
        System.out.println(message);
    }
}

