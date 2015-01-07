package com.Sts.UI;

import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/** This is a dialog centered on a parent with a displayPanel displayed and an okCancelObject which it communicates with
 *  when ok/apply/cancel buttons are pushed.
 */


public class StsOkCancelDialog extends JDialog
{
	public StsDialogFace[] okCancelObjects;
	private StsGroupBox buttonBox = new StsGroupBox();
	private StsButton okButton = new StsButton("OK", "Accept changes and dismiss dialog.", this, "ok");
	private StsButton cancelButton = new StsButton("Cancel", "Reject changes and dismiss dialog.", this, "cancel");

	public StsOkCancelDialog(Frame frame, StsDialogFace[] okCancelObjects, String title, boolean modal)
	{
		super(frame, title, modal);
        if(frame != null) setLocation(frame.getLocation());
		this.okCancelObjects = okCancelObjects;
		layoutPanels();
		addWindowCloseOperation();
		display();
	}

	public StsOkCancelDialog(Frame frame, StsDialogFace okCancelObject, String title, boolean modal)
	{
		this(frame, new StsDialogFace[] { okCancelObject }, title, modal);
	}

	public void display()
	{
		pack();
		setVisible(true);
	}

	private void layoutPanels()
	{
		StsJPanel panel = StsJPanel.addInsets();
		panel.gbc.fill = GridBagConstraints.HORIZONTAL;
		for(int n = 0; n < okCancelObjects.length; n++)
			panel.add(okCancelObjects[n].getPanel());
		Insets insets = new Insets(4, 4, 4, 4);
		okButton.setMargin(insets);
		cancelButton.setMargin(insets);
		buttonBox.addToRow(okButton);
		buttonBox.addEndRow(cancelButton);
		panel.add(buttonBox);
		getContentPane().add(panel);
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
					 close();
				 }
			 }
		);
	}

	public void close()
	{
		for(int n = 0; n < okCancelObjects.length; n++)
			okCancelObjects[n].dialogSelectionType(StsDialogFace.CLOSING);
		dispose();
	}

	public void ok()
	{
		for(int n = 0; n < okCancelObjects.length; n++)
			okCancelObjects[n].dialogSelectionType(StsDialogFace.OK);
		dispose();
        if(okCancelObjects.length == 1) setVisible(false);
    }

	public void cancel()
	{
		for(int n = 0; n < okCancelObjects.length; n++)
			okCancelObjects[n].dialogSelectionType(StsDialogFace.CANCEL);
		setVisible(false);
		dispose();
	}

	static public void main(String[] args)
	{
		StsDialogFaceStringPanel testObject = new StsDialogFaceStringPanel("Test OK-Cancel String", "String", "string");
		StsOkCancelDialog dialog = new StsOkCancelDialog((Frame)null, new StsDialogFace[] { testObject }, "Title.", true);
	}
}

