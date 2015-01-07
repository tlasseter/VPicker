
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

public class StsPasswordDialog extends JDialog
{
    private String text = null;

	private JPanel mainPanel = new JPanel();
	private BorderLayout paneLayout = new BorderLayout();
	private JPanel buttonPanel = new JPanel();
	protected JButton okayButton = new JButton();
	private FlowLayout flowLayout1 = new FlowLayout();
	private TextField textField = new TextField();

	public StsPasswordDialog(Frame parent, String title, String text)
	{
        super(parent, title, true);
        this.setLocationRelativeTo(parent);
        this.text = text;
		jbInit();
		pack();
	}

	private void jbInit()
	{
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		mainPanel.setLayout(paneLayout);
		buttonPanel.setLayout(flowLayout1);
		okayButton.setText("Okay");

		okayButton.addActionListener(new java.awt.event.ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				okayButton_actionPerformed(e);
			}
		});

		buttonPanel.add(okayButton);
		mainPanel.add(buttonPanel, BorderLayout.SOUTH);

        textField.setText(text);
        textField.setEchoChar('*');

		mainPanel.add(textField, BorderLayout.NORTH);
        getContentPane().add(mainPanel);
	}

	protected void okayButton_actionPerformed(ActionEvent e)
	{
		dispose();
    }

    /** get/set the text */
    public void setText(String text) { textField.setText(text); }
    public String getText() { return textField.getText(); }

    public static void main(String[] args)
    {
        StsPasswordDialog d = new StsPasswordDialog(null, "StsTextAreaDialog test",
                "Display text:\n2nd line\n3rd line");
        d.setVisible(true);
        System.out.println("StsPasswordDialog:  text = " + d.getText());
    }

}

