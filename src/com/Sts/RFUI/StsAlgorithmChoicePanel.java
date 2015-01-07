
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.RFUI;

import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;

public class StsAlgorithmChoicePanel extends JPanel
{
    StsGroupBox groupBox1 = new StsGroupBox();
    GridLayout gridLayout1 = new GridLayout();
    GridLayout gridLayout2 = new GridLayout();
	ButtonGroup group = new ButtonGroup();
    JRadioButton jRadioButton1 = new JRadioButton();
    JRadioButton jRadioButton2 = new JRadioButton();
    JRadioButton jRadioButton3 = new JRadioButton();
    private boolean invokedStandalone = false;

    public StsAlgorithmChoicePanel()
    {
        try
        {
            jbInit();
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

    private void jbInit() throws Exception
    {
        group.add(jRadioButton1);
        group.add(jRadioButton2);
        group.add(jRadioButton3);

        gridLayout2.setColumns(1);
        gridLayout2.setRows(0);

        jRadioButton1.setSelected(true);
        jRadioButton1.setText("Inverse Distance");
        jRadioButton2.setText("Krige");
        jRadioButton3.setText("Sequential Gaussian Simulation");

        setLayout(gridLayout1);
        add(groupBox1, null);

        groupBox1.setLabel("Algorithm Type");
        groupBox1.setLayout(gridLayout2);
        groupBox1.add(jRadioButton1, null);
//        groupBox1.add(jRadioButton2, null);
        groupBox1.add(jRadioButton3, null);
    }

    public int getChoice()
    {
        if( jRadioButton1.isSelected() ) return 1;
        if( jRadioButton2.isSelected() ) return 2;
        if( jRadioButton3.isSelected() ) return 3;
        return 0;
    }

    public static void main(String[] args)
    {
        StsAlgorithmChoicePanel StsAlgorithmChoicePanel = new StsAlgorithmChoicePanel();
        StsAlgorithmChoicePanel.invokedStandalone = true;

        JFrame frame = new JFrame();
        frame.getContentPane().add(StsAlgorithmChoicePanel);
        frame.pack();
        frame.setVisible(true);
    }

}
