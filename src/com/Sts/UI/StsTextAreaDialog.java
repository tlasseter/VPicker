
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// dialog box to set a text field

package com.Sts.UI;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class StsTextAreaDialog extends JDialog
{
    private String text = null;
    private int rows;
    private int columns;
    private boolean scrollbarsOn;

	private JPanel mainPanel = new JPanel();
	private BorderLayout paneLayout = new BorderLayout();
	private JPanel buttonPanel = new JPanel();
	protected JButton okayButton = new JButton();
	//protected JButton cancelButton = new JButton();
	//protected JButton helpButton = new JButton();
	private FlowLayout flowLayout1 = new FlowLayout();
	private JTextArea textArea = new JTextArea();
	private JScrollPane textPane;

	public StsTextAreaDialog(Frame parent, String title, String text)
	{
        this(parent, title, text, 20, 30);
    }

	public StsTextAreaDialog(Frame parent, String title, boolean modal)
	{
		this(parent, title, null, 20, 30, modal);
    }

	public StsTextAreaDialog(Frame parent, String title, String text, int rows, int columns)
	{
        this(parent, title, text, rows, columns, true);
    }

	public StsTextAreaDialog(Frame parent, String title, String text, int rows, int columns, boolean modal)
	{
        this(parent, title, text, rows, columns, modal, false);
    }

	public StsTextAreaDialog(Frame parent, String title, String text, int rows, int columns, boolean modal, boolean scrollbarsOn)
	{
        super(parent, title, modal);
        this.setLocationRelativeTo(parent);
        this.rows = rows;
        this.text = text;
        this.columns = columns;
        this.scrollbarsOn = scrollbarsOn;
		jbInit();
		pack();
	}

	private void jbInit()
	{
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		mainPanel.setLayout(paneLayout);
		buttonPanel.setLayout(flowLayout1);
		okayButton.setText("Okay");
		//cancelButton.setText("Cancel");
		//helpButton.setText("Help");

		okayButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				okayButton_actionPerformed(e);
			}
		});
/*
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
*/
		buttonPanel.add(okayButton);
		//buttonPanel.add(cancelButton);
		//buttonPanel.add(helpButton);
		mainPanel.add(buttonPanel, BorderLayout.SOUTH);

        textArea.setRows(rows);
        textArea.setColumns(columns);
        textArea.setText(text);

        if (!scrollbarsOn) textPane = new JScrollPane(textArea);
        else textPane = new JScrollPane(textArea, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

		mainPanel.add(textPane, BorderLayout.NORTH);
        getContentPane().add(mainPanel);
	}

	protected void okayButton_actionPerformed(ActionEvent e)
	{
		dispose();
   }

/*
	protected void cancelButton_actionPerformed(ActionEvent e)
	{
		dispose();
	}

	protected void helpButton_actionPerformed(ActionEvent e)
	{

	}
*/

    /** get/set the text */
    public void setText(String text) { textArea.setText(text); }
    public String getText() { return textArea.getText(); }

    public void appendLine(String line)
    {
       append(line + "\n");
    }

    public void append(String string)
    {
        textArea.append(string);
    }

    public static void main(String[] args)
    {
        StsTextAreaDialog d = new StsTextAreaDialog(null, "StsTextAreaDialog test",
                "Display text:\n2nd line\n3rd line");
        d.setVisible(true);
        System.out.println("StsTextAreaDialog:  text = " + d.getText());
    }

}

