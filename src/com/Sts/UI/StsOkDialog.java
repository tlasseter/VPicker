package com.Sts.UI;

import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;


/** This is a dialog centered on a parent with a panel added and OK button below. */
public class StsOkDialog extends JDialog
{
    StsJPanel backPanel = StsJPanel.addInsets();
	private StsGroupBox buttonBox = new StsGroupBox();
	private StsButton okButton = new StsButton("OK", "Accept Changes and Exit", this, "ok");
    transient boolean okPressed = false;

    public StsOkDialog() {}
    public StsOkDialog(Frame frame, JPanel panel, String title, boolean modal)
	{
		super(frame, title, modal);
        if(frame != null) StsToolkit.centerComponentOnFrame(this, frame);
        constructPanel(panel);
		addWindowCloseOperation();
		display();
	}

    public void display()
	{
		pack();
		setVisible(true);
	}

    private void constructPanel(JPanel panel)
	{
		backPanel.add(panel);
		Insets insets = new Insets(4, 4, 4, 4);
		okButton.setMargin(insets);
		buttonBox.addToRow(okButton);
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
        okPressed = false;
		dispose();
	}

	public void ok()
	{
        okPressed = true;
        dispose();
	}

    public boolean okPressed() { return okPressed; }
}