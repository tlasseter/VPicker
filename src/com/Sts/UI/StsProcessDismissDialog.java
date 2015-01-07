package com.Sts.UI;

import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */

/** This is a dialog centered on a parent with a displayPanel displayed and an okCancelObject which it communicates with
 *  when ok/apply/cancel buttons are pushed.
 */


public class StsProcessDismissDialog extends JDialog
{
	private StsDialogFace processObject;
	private StsGroupBox buttonBox = new StsGroupBox();
	private StsButton processButton = new StsButton("Process", "Execute this process.", this, "process");
	private StsButton cancelButton = new StsButton("Cancel", "Quit without processing.", this, "cancel");

	public StsProcessDismissDialog(Frame frame, Component displayPanel, StsDialogFace processObject, String title, boolean modal)
	{
		super(frame, title, modal);
		if (frame != null)
			setLocation(frame.getLocation());
		layoutPanel(displayPanel);
		this.processObject = processObject;
		addWindowCloseOperation();
		pack();
		setVisible(true);
	}

	public StsProcessDismissDialog(Frame frame, StsDialogFace processObject, String title, boolean modal)
	{
		super((Frame)null, title, modal);
		layoutPanel(processObject.getPanel());
		this.processObject = processObject;
		addWindowCloseOperation();
		pack();
		setVisible(true);
	}

	private void layoutPanel(Component displayPanel)
	{
		StsJPanel panel = StsJPanel.addInsets();
		panel.gbc.fill = GridBagConstraints.HORIZONTAL;
		if(displayPanel != null) panel.add(displayPanel);
		Insets insets = new Insets(4, 4, 4, 4);
		processButton.setMargin(insets);
		cancelButton.setMargin(insets);
		buttonBox.addToRow(processButton);
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
					 dispose();
					 processObject.dialogSelectionType(StsDialogFace.CLOSING);
				 }
			 }
		);
	}

	public void process()
	{
		processObject.dialogSelectionType(StsDialogFace.PROCESS);
		dispose();
	}

	public void cancel()
	{
		processObject.dialogSelectionType(StsDialogFace.CANCEL);
		dispose();
	}

	static public void main(String[] args)
	{
		TestDismissObject testObject = new TestDismissObject();
		StsProcessDismissDialog dialog = new StsProcessDismissDialog((Frame)null, testObject, "Title.", true);
	}
}

class TestProcessDismiss implements StsDialogFace
{
	String string = "test";
	StsStringFieldBean stringBean = new StsStringFieldBean(this, "string", "String");

	TestProcessDismiss()
	{
	}

    public StsDialogFace getEditableCopy()
    {
        return (StsDialogFace) StsToolkit.copyObjectNonTransientFields(this);
    }

    public void dialogSelectionType(int type)
	{
		System.out.println("Selection Type " + type);
	}
	public Component getPanel(boolean val) { return getPanel(); }
	public Component getPanel()
	{
		StsGroupBox groupBox = new StsGroupBox("Test Process/Dismiss");
		groupBox.add(stringBean);
		return (Component)groupBox;
	}

	public void setString(String s) { string = s; }
	public String getString() { return string; }
}

class TestDismissObject implements StsDialogFace
{
	String string = "test";
	StsStringFieldBean stringBean = new StsStringFieldBean(this, "string", "String");
	TestDismissObject()
	{
	}

	public void dialogSelectionType(int type)
	{
		System.out.println("Selection Type " + type);
	}
	public Component getPanel(boolean val) { return getPanel(); }
	public Component getPanel()
	{
		StsGroupBox groupBox = new StsGroupBox("Test OkCancel");
		groupBox.add(stringBean);
		return (Component)groupBox;
	}

    public StsDialogFace getEditableCopy()
    {
        return (StsDialogFace) StsToolkit.copyObjectNonTransientFields(this);
    }

    public void setString(String s) { string = s; }
	public String getString() { return string; }
}
