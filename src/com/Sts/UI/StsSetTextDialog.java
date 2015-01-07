
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// dialog box to set a text field

package com.Sts.UI;

import com.Sts.Help.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class StsSetTextDialog extends JDialog
{
    static final private String HELP_URL = "settextdialog.html";
    static final private String HELP_TARGET = "ui.setTextDialog";

    private String text = null;
    private int columns;
    private String relativeHelpFile = HELP_URL;

    private JFrame parent;
	private JPanel mainPanel = new JPanel();
	private BorderLayout paneLayout = new BorderLayout();
	private JPanel buttonPanel = new JPanel();
	protected JButton okayButton = new JButton();
	protected JButton cancelButton = new JButton();
	protected JButton helpButton = new JButton();
	private FlowLayout flowLayout1 = new FlowLayout();
	protected JTextField textField = new JTextField();

	public StsSetTextDialog(JFrame parent, String title, String initialText)
	{
        this(parent, title, initialText, 25);
    }
	public StsSetTextDialog(JFrame parent, String title, String initialText,
            int columns)
	{
        super(parent, title, true);  // modal
        this.parent = parent;
        setText(initialText);
        this.columns = columns;
		jbInit();
		pack();
	}

	private void jbInit()
	{
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		mainPanel.setLayout(paneLayout);
		buttonPanel.setLayout(flowLayout1);
		okayButton.setText("Okay");
		cancelButton.setText("Cancel");
		helpButton.setText("Help");

		okayButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				okayButton_actionPerformed(e);
			}
		});
		cancelButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				cancelButton_actionPerformed(e);
			}
		});
		helpButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				helpButton_actionPerformed(e);
			}
		});
        setHelpTarget(HELP_TARGET);

		buttonPanel.add(okayButton);
		buttonPanel.add(cancelButton);
		buttonPanel.add(helpButton);
		mainPanel.add(buttonPanel, BorderLayout.SOUTH);

        textField.setColumns(columns);
		mainPanel.add(textField, BorderLayout.NORTH);
        getContentPane().add(mainPanel);
	}

	protected void okayButton_actionPerformed(ActionEvent e)
	{
        text = textField.getText();
		dispose();
   }

	protected void extraButton_actionPerformed(ActionEvent e)
	{
        // override in subclass
    }

	protected void cancelButton_actionPerformed(ActionEvent e)
	{
        text = null;
		dispose();
	}

	protected void helpButton_actionPerformed(ActionEvent e)
	{
//        HelpManager.loadRelativeURL(relativeHelpFile);
	}

    /** set/get the text */
    public void setText(String text)
    {
        textField.setText(text);
        this.text = text;
    }
    public String getText()
    {
        return text;
    }

    /** set the help link (once) */
    public boolean setHelpTarget(String target)
    {
        return HelpManager.setButtonHelp(helpButton, HelpManager.GENERAL, target, parent);
    }
    public void setHelpFile(String relativeFile) { relativeHelpFile = relativeFile; }


    public static void main(String[] args)
    {
        StsSetTextDialog d = new StsSetTextDialog(null, "StsSetTextDialog test",
                "Change me");
        d.setVisible(true);
        System.out.println("StsSetTextDialog:  text set = " + d.getText());
        System.exit(0);
    }

}

