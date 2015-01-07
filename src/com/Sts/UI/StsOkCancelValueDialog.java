package com.Sts.UI;

import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/** This is a dialog centered on a parent with a displayPanel displayed and an okCancelObject which it communicates with
 *  when ok/apply/cancel buttons are pushed.
 */


public class StsOkCancelValueDialog extends JDialog
{
    StsJPanel mainPanel;
    private float value;
    private StsFloatFieldBean valueBean;
	private StsGroupBox buttonBox = new StsGroupBox();
	private StsButton okButton = new StsButton("OK", "Accept changes and dismiss dialog.", this, "ok");
	private StsButton cancelButton = new StsButton("Cancel", "Reject changes and dismiss dialog.", this, "cancel");

	public StsOkCancelValueDialog(Frame frame, String title, String valueLabel)
	{
		super(frame, title, true);
        if(frame != null) setLocation(frame.getLocation());
		constructMainPanel();
        addValueBean(valueLabel);
        addWindowCloseOperation();
		display();
	}

	public StsOkCancelValueDialog(Frame frame, String title, StsJPanel panel, String valueLabel)
	{
		super(frame, title, true);
        if(frame != null) setLocation(frame.getLocation());
        constructMainPanel();
        addPanel(panel);
        addValueBean(valueLabel);
		addWindowCloseOperation();
		display();
	}

    public void display()
	{
		pack();
		setVisible(true);
	}

    private void constructMainPanel()
    {
        mainPanel = StsJPanel.addInsets();
        mainPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
    }

    private void addPanel(StsJPanel panel)
    {
        mainPanel.add(panel);
    }

    private void addValueBean(String valueLabel)
    {
        valueBean = new StsFloatFieldBean(this, "value", valueLabel);
        mainPanel.add(valueBean);
		Insets insets = new Insets(4, 4, 4, 4);
		okButton.setMargin(insets);
		cancelButton.setMargin(insets);
		buttonBox.addToRow(okButton);
		buttonBox.addEndRow(cancelButton);
		mainPanel.add(buttonBox);
		getContentPane().add(mainPanel);
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
		dispose();
	}

	public void ok()
	{
		dispose();
    }

	public void cancel()
	{
		setVisible(false);
		dispose();
	}

	static public void main(String[] args)
	{
		StsOkCancelValueDialog dialog = new StsOkCancelValueDialog((Frame)null, "Title.", "Change value");
	}

    public float getValue()
    {
        return value;
    }

    public void setValue(float value)
    {
        this.value = value;
    }
}
