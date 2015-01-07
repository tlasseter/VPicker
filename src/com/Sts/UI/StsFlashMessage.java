package com.Sts.UI;

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Sep 5, 2007
 * Time: 10:02:37 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsFlashMessage implements Runnable
{
    Frame owner;
    String title;
    String message;
    int delaySecs;
    JDialog dialog;

    public StsFlashMessage(Frame owner, String message)
    {
        this(owner, "Flash message", message, 2);
    }

    public StsFlashMessage(Frame owner, String message, int delaySecs)
    {
        this(owner, "Flash message", message, delaySecs);
    }

    public StsFlashMessage(Frame owner_, String title_, String message_, int delaySecs_)
    {
        this.owner = owner_;
        this.title = title_;
        this.message = message_;
        this.delaySecs = delaySecs_;
        StsToolkit.runWaitOnEventThread(this);
    }

    public void run()
    {
        dialog = new JDialog(owner, title, false);
        JTextField textField = new JTextField(message);
        dialog.getContentPane().add(textField);
        dialog.pack();
        if(owner != null) StsToolkit.centerComponentOnFrame(dialog, owner);
        dialog.setVisible(true);
    }

    static public void constructor(Frame frame, String message)
    {
        StsFlashMessage flashMessage = new StsFlashMessage(frame, message);

        try { Thread.currentThread().sleep(2000); }
		catch(Exception e) { }

        flashMessage.dispose();
    }

    public void dispose()
    {
        if(dialog != null) dialog.setVisible(false);
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
        StsFlashMessage.constructor(frame, "Flash message");
/*
        StsFlashMessage flashMessage = new StsFlashMessage(frame, "Flash message");

        try { Thread.currentThread().sleep(2000); }
		catch(Exception e) { }

        flashMessage.dispose();
*/
        frame.dispose();
    }
}
