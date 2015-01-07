package com.Sts.UI;

import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version 1.0
 */

public class StsProgressDialog extends JDialog
{
	private StsProgressPanel progressPanel = null;

	public StsProgressDialog(JFrame frame, String title, boolean modal)
	{
		super(frame, title, modal);
		if (frame != null)
        {
            setLocation(frame.getLocation());
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        }
		progressPanel = StsProgressPanel.constructorWithCancelButton(5, 50);
		layoutPanel();
		addWindowCloseOperation();
		pack();
		setVisible(true);
	}

	public StsProgressPanel getProgressPanel()
	{
		return progressPanel;
	}

	private void layoutPanel()
	{
		StsJPanel panel = StsJPanel.addInsets();
		panel.gbc.fill = GridBagConstraints.BOTH;
		panel.gbc.weighty = 1;
		panel.addEndRow(progressPanel);
		panel.gbc.anchor = GridBagConstraints.SOUTH;
		getContentPane().add(panel);
		setResizable(false);
	}

	private void addWindowCloseOperation()
	{
		setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
		addWindowListener
		 (
			 new WindowAdapter()
			 {
				 public void windowClosing(WindowEvent e)
				 {
					 progressPanel.cancel();
                     dispose();
                 }
			 }
		);
	}
}
