
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;


import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class StsRadioButtonSelectPanel extends JPanel
{
	GridBagLayout gridBagLayout = new GridBagLayout();
	JLabel title = new JLabel();
    ButtonGroup group = new ButtonGroup();
	JRadioButton[] buttons;
	int nButtons;

	static final int nTitleRows = 3;

	public StsRadioButtonSelectPanel()
    {
    }

    public StsRadioButtonSelectPanel(String title, String[] buttonLabels)
    {
    	setTitle(title);
		nButtons = buttonLabels.length;
		buttons = new JRadioButton[nButtons];
		for(int n = 0; n < nButtons; n++)
			buttons[n] = new JRadioButton(buttonLabels[n]);
		try { jbInit(); }
        catch(Exception e) { e.printStackTrace(); }
    }


    private void jbInit()
    {
		this.setLayout(gridBagLayout);

		title.setHorizontalAlignment(0);
		int row = 0;
		this.add(title, new GridBagConstraints(0, 0, 1, nTitleRows, 1.0, 0.2
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		row += nTitleRows;
		for(int n = 0; n < nButtons; n++)
		{
			JRadioButton button = buttons[n];
		    button.setHorizontalAlignment(0);

			this.add(button, new GridBagConstraints(0, row++, 1, 1, 1.0, 0.0
				,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 10, 0, 0), 0, 0));

            group.add(button);
		}
		buttons[0].setSelected(true);
    }


	public void setButtonIndexSelected(int buttonIndexSelected)
	{
        if(buttons == null || buttons.length <= buttonIndexSelected) return;
        buttons[buttonIndexSelected].setSelected(true);
	}

	public int getButtonIndexSelected()
	{
		for(int n = 0; n < nButtons; n++)
		    if(buttons[n].isSelected()) return n;
		return 0;
	}

	public boolean isButtonSelected(int index)
	{
		return buttons[index].isSelected();
	}

    public void setTitle(String s) { title.setText(s); }

    public void addButtonListeners(ItemListener listener)
    {
        if(buttons == null) return;
        int nButtons = buttons.length;
        for(int n = 0; n < nButtons; n++)
            buttons[n].addItemListener(listener);
    }
}
