package com.Sts.MVC.ViewWell;

/**
 * Title:        Sts Well Viewer Description:  Well Model-Viewer Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class StsWellWindowFrame extends StsWindow
{
    String title;
    StsWellFrameViewModel wellViewModel;
    JMenuBar operationsMenuBar = new JMenuBar();
    JMenu operationsMenu = new JMenu("Select Operation");

    public StsWellWindowFrame(StsWellFrameViewModel wellViewModel, Point location)
    {
        setLocation(location);
        this.wellViewModel = wellViewModel;
        setTitle(wellViewModel.getName());
        addWindowListener();
        actionManager = new StsActionManager(model); // jbw need model
    }

    public void start()
    {
        validate();
        pack();
        wellViewModel.initAndStartGL();
        setVisible(true);
    }

    //	public StsWellTrackGLPanel getWellTrackDisplayPanel() { return wellTrackDisplayPanel; }
    public void rebuild()
    {
//		setSize(new Dimension(getWidth(), getHeight()));
        validate();
        pack();
        setVisible(true);
    }

    public void printWindow()
    {
    }

    private void addWindowListener()
    {
        super.addWindowListener(
            new WindowAdapter()
            {
                public void windowClosing(WindowEvent e)
                {
                    wellViewModel.closeWindow();
                }
                public void windowStateChanged(WindowEvent e)
                 {
                     System.out.println(e.toString());
                 }
            });
		super.addComponentListener(
			  new ComponentAdapter()
		     {
			    public void componentMoved(ComponentEvent e)
				{
					Component c = e.getComponent();
					wellViewModel.setFrameLocation(c.getLocation());
				}
				public void componentResized(ComponentEvent e)
				{
					Component c = e.getComponent();
					wellViewModel.setFrameSize(c.getSize());
				}

		     });
    }
/*
	public void paint(Graphics g)
	{
		super.paint(g);
		if (wellViewModel != null)
		{
			wellViewModel.setFrameSize(getSize());
			wellViewModel.setFrameLocation(getLocation());
		}
	}
*/
}
