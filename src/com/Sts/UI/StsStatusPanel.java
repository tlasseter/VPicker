
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;

public class StsStatusPanel extends JPanel implements StsStatusUI
{
    static public int UI_UPDATE_DELAY = 5;
    Border border = BorderFactory.createLoweredBevelBorder();
    GridBagLayout gridBag = new GridBagLayout();
    JProgressBar progress = new JProgressBar();
    JLabel titleBar = new JLabel();
    JLabel status = new JLabel();

    static JDialog dialog = null;

    public StsStatusPanel()
    {
        super(false);

        try { jbInit(); }
        catch(Exception e) { }
    }

    private void jbInit() throws Exception
    {
        setLayout(gridBag);
        Color bg = progress.getBackground();
        setBackground(bg);
		status.setBackground(bg);
        status.setMaximumSize(new Dimension(400, 25));
        status.setMinimumSize(new Dimension(400, 25));
        status.setPreferredSize(new Dimension(400, 25));

		titleBar.setBackground(bg);
//		logMessage("status");
		status.setHorizontalAlignment(0);

		this.add(titleBar, new GridBagConstraints(0, 0, 1, 1, 0.0, 1.0
            ,GridBagConstraints.NORTH, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 0, 0), 0, 0));
		this.add(status, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 20, 5, 20), 0, 0));
		this.add(progress, new GridBagConstraints(0, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(20, 20, 20, 20), 0, 0));
	}

    public void setMaximum(float max) { setMaximum((int)max); }
    public void setMinimum(float min) { setMinimum((int)min); }
    public void setProgress(float n) { setProgress((int)n); }

    public void setMaximum(int max) { progress.setMaximum(max); }
    public void setMinimum(int min) { progress.setMinimum(min); }
    public void setProgress(int n_)
    {
        final int n = n_;
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    progress.setValue(n);
                    progress.repaint();
                }
            }
        );
    }

    public int getProgress()
    {
    	return progress.getValue();
    }
    public void setTitle(String msg_)
    {
        final String msg = msg_;
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    titleBar.setText(msg);
                    titleBar.repaint();
                }
            }
        );
    }

	public void incrementProgress()
	{
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    int progressValue = progress.getValue();
                    progress.setValue(++progressValue);
                }
            }
        );
	}

    public void setText(String msg) { setText(msg, 0, true); }
    public void setText(String msg, int msec){ setText(msg, msec, true); }
    public void setText(String msg, boolean log) { setText(msg, 0, log); }
    public void setText(String msg_, int msec, boolean log_)
    {
        final String msg = msg_;       
        if(log_) StsMessageFiles.logMessage(msg);
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    status.setText(msg);
                }
            }
        );
    }

    public void sleep(int msec)
    {
        if (msec<1) msec = 1;
        if (msec>10000) msec = 10000;
        try { Thread.currentThread().sleep(msec); }
        catch(Exception e) { }
    }

    static public StsStatusPanel constructStatusDialog(Frame frame, String title)
    {

		dialog = new JDialog(frame, title, false);
        StsStatusPanel status = new StsStatusPanel();
        dialog.getContentPane().add(status);
        dialog.setLocationRelativeTo(frame);
        dialog.pack();
        dialog.setVisible(true);

        return status;
    }

    static public void disposeDialog()
    {
        if(dialog == null) return;
        dialog.dispose();
    }

    static public void main(String[] args)
    {
        StsStatusPanel statusPanel = constructStatusDialog((Frame)null, "Status Panel");

        statusPanel.setTitle("title of panel .....");
        statusPanel.setMaximum(1000);
        statusPanel.setMinimum(0);

        for( int i=0; i<1000; i++ )
        {
            StsMessageFiles.logMessage("Working on line " + i);
            statusPanel.setProgress(i);
            try { Thread.currentThread().sleep(1); }
            catch(Exception e) { }
        }
    }
}
