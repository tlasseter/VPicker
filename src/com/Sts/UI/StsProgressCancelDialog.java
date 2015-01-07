package com.Sts.UI;

import com.Sts.Interfaces.*;
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

public class StsProgressCancelDialog extends JDialog
{
	private boolean processCompleted = false;
	private StsDialogFace processObject;
	private StsButton cancelButton = new StsButton("Cancel", "Quit processing.", this, "cancel");
	private StsProgressPanel progressPanel;

	public StsProgressCancelDialog(JFrame frame, StsDialogFace processObject, String title, String progressBarTitle,boolean modal)
	{
		super(frame, title, modal);
		if (frame != null)
        {
            setLocation(frame.getLocation());
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        }
        progressPanel = StsProgressPanel.constructor(10, 40);
		layoutPanel();
		this.processObject = processObject;
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
		JPanel buttonPanel = new JPanel();
		panel.gbc.weighty = 0;
		buttonPanel.add(cancelButton);
		panel.addEndRow(buttonPanel);
		getContentPane().add(panel);
		setResizable(true);
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
					 cancel();
				 }
			 }
		);
	}

	public void cancel()
	{
		if (processCompleted)
		{
			processObject.dialogSelectionType(StsDialogFace.DISMISS);
			cancelButton.setEnabled(false);
			dispose();
        }
        else
        {
            int answer = JOptionPane.showConfirmDialog(this, "Are you sure you want to cancel the process?", "Cancel Process", JOptionPane.YES_NO_OPTION);
            if (answer == JOptionPane.YES_OPTION)
            {
                processObject.dialogSelectionType(StsDialogFace.CANCEL);
                cancelButton.setEnabled(false);
                dispose();
            }
        }
    }

	public void setProcessCompleted()
	{
		cancelButton.setText("Close");
		cancelButton.setToolTipText("Dismiss this dialog");
		processCompleted = true;
	}
}
