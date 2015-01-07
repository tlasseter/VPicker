

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;

public class StsSectionGapDialog extends JDialog
	implements ActionListener
{
	static final String TITLE = "Section Gaps";
	static final String OK_LABEL = "OK";
	static final String CANCEL_LABEL = "Cancel";
	static final String OK_COMMAND = "OK";
	static final String CANCEL_COMMAND = "Cancel";

    Border lineBorder = LineBorder.createBlackLineBorder();
    BevelBorder border = (BevelBorder)BorderFactory.createLoweredBevelBorder();
	GridBagLayout gridBag = new GridBagLayout();
	FlowLayout flow = new FlowLayout();
    JPanel contentPane = new JPanel();
    JPanel topPanel = new JPanel();
	JPanel buttonPanel = new JPanel();
    JPanel solidBox = new JPanel();
	JPanel lightBox = new JPanel();
	JTextField lightTextField = new JTextField();
	JTextField solidTextField = new JTextField();
	JButton okButton = new JButton();
	JButton cancelButton = new JButton();
	BorderLayout borderLayout = new BorderLayout();
    private String buttonPressed = null;

	public StsSectionGapDialog()
	{
		try
		{
			jbInit();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	private void jbInit() throws Exception
	{
		solidBox.setPreferredSize(new Dimension(20, 20));
		solidBox.setBackground(new Color(0, 0, 95));
		solidBox.setBorder(lineBorder);
		lightBox.setPreferredSize(new Dimension(20, 20));
		lightBox.setBackground(Color.blue);
		lightBox.setBorder(lineBorder);

		okButton.setText(OK_LABEL);
		okButton.setActionCommand(OK_COMMAND);
		cancelButton.setText(CANCEL_LABEL);
		cancelButton.setActionCommand(CANCEL_COMMAND);

		okButton.addActionListener(this);
		cancelButton.addActionListener(this);


		topPanel.setLayout(gridBag);
		topPanel.setBorder(border);
		topPanel.add(lightBox, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
		topPanel.add(solidBox, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
		topPanel.add(lightTextField, new GridBagConstraints(1, 0, 2, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 5), 0, 0));
		topPanel.add(solidTextField, new GridBagConstraints(1, 1, 2, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 5), 0, 0));

		flow.setAlignment(FlowLayout.RIGHT);
		buttonPanel.setLayout(flow);
		buttonPanel.add(okButton);
		buttonPanel.add(cancelButton);

		contentPane.setLayout(borderLayout);
		contentPane.setPreferredSize(new Dimension(200, 100));
		setContentPane(contentPane);
        contentPane.add(topPanel, BorderLayout.CENTER);
        contentPane.add(buttonPanel, BorderLayout.SOUTH);


        setModal(true);
		setResizable(false);
		setTitle(TITLE);
        pack();
	}

	public static void main(String[] args)
	{
		StsSectionGapDialog d = new StsSectionGapDialog();
        d.setRightGap(2.5f);
        d.setLeftGap(1.0f);
        d.setVisible(true);

        System.out.println("Gap values> front: " + d.getRightGap() + " | back: " + d.getLeftGap());
        d.dispose();
        System.exit(0);
	}

	public void setGapColor(Color color)
    {
    	solidBox.setBackground(color.darker());
    	lightBox.setBackground(color);
    }

    public void setRightGap(float gap)
    {
		solidTextField.setText(String.valueOf(gap));
    }

    public float getRightGap()
    {
		String gapStr = solidTextField.getText();
        Float f = new Float(gapStr);
        return f.floatValue();
    }
    public void setLeftGap(float gap)
    {
		lightTextField.setText(String.valueOf(gap));
    }

    public float getLeftGap()
    {
		String gapStr = lightTextField.getText();
        Float f = new Float(gapStr);
        return f.floatValue();
    }

	public void actionPerformed(ActionEvent e)
	{
    	if( e.getSource() instanceof JButton )
        {
        	setVisible(false);
            buttonPressed = e.getActionCommand();
        }
	}

    public boolean okWasPressed()
    {
		if( buttonPressed == null )
        	return false;
        else
        	return buttonPressed.equals(OK_COMMAND);
    }


}


