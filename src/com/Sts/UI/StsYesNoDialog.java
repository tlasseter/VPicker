
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/** Use JOptionPane to display a message of type WARNING, INFO, or FLASH.  A
 *  FLASH message stays up for a few seconds and automatically disappears.  The
 *  other two types have an OK button to dismiss.  The message shouldn't have
 *  more than 32 characters or so per line or it will be truncated by JOptionPane:
 *  Insert a "\n" to make a new line.
 */
public class StsYesNoDialog implements Runnable
{
	private Component component;
	private String messageString;
    public boolean answer;


	static public boolean questionValue(Component component, String string)
	{
        StsYesNoDialog yesNoQuestion = new StsYesNoDialog(component, string);
        return yesNoQuestion.answer;
	}

	public StsYesNoDialog(Component component, String messageString)
	{
		this.component = component;
		this.messageString = messageString;
        StsToolkit.runWaitOnEventThread(this);
	}

    public void run()
	{
        int yesOrNo = JOptionPane.showConfirmDialog(component, messageString,
                        "Question Message", JOptionPane.YES_NO_OPTION,
                        JOptionPane.QUESTION_MESSAGE);

        answer = yesOrNo == JOptionPane.YES_OPTION;
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

        String message = "Do you want to do this?";
        boolean answer = questionValue(frame, message);
    }
}