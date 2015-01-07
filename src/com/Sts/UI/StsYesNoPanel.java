
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;


import javax.swing.*;
import java.awt.*;

public class StsYesNoPanel extends JPanel
{
	GridBagLayout gridBagLayout = new GridBagLayout();
	JLabel questionLabel = new JLabel();
    ButtonGroup group = new ButtonGroup();
	JRadioButton yesButton = new JRadioButton();
	JRadioButton noButton = new JRadioButton();
    Object[] items = null;

	static final int nTitleRows = 3;

	public StsYesNoPanel()
    {
    	this("This is the question?", "Yes", "No");
    }

    public StsYesNoPanel(String question, String yesText, String noText)
    {
    	setQuestion(question);
        setYesText(yesText);
        setNoText(noText);
		try { jbInit(); }
        catch(Exception e) { e.printStackTrace(); }
    }


    private void jbInit()
    {
		this.setLayout(gridBagLayout);
		questionLabel.setHorizontalAlignment(0);
		yesButton.setHorizontalAlignment(0);
		noButton.setHorizontalAlignment(0);
		int rows = 0;
		this.add(questionLabel, new GridBagConstraints(0, 0, 1, nTitleRows, 1.0, 0.2
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		rows += nTitleRows;
		this.add(yesButton, new GridBagConstraints(0, rows, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 10, 0, 0), 0, 0));
		rows += 1;
		this.add(noButton, new GridBagConstraints(0, rows, 1, 1, 1.0, 1.0
            ,GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 10, 0, 0), 0, 0));
        yesButton.setSelected(true);
        group.add(yesButton);
        group.add(noButton);
    }

    public void setQuestion(String s) { questionLabel.setText(s); }
    public String getQuestion() { return questionLabel.getText(); }
    public void setYesText(String s) { yesButton.setText(s); }
    public String getYesText() { return yesButton.getText(); }

    public boolean isYes() { return yesButton.isSelected(); }
    public boolean isNo() { return noButton.isSelected(); }

    public void setNoText(String s) { noButton.setText(s); }
    public String getNoText() { return noButton.getText(); }
}
