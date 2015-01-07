
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
public class StsPowerPanel extends JPanel
{
    GridLayout gridLayout1 = new GridLayout();
    JLabel jLabel1 = new JLabel();
    StsGroupBox groupBox1 = new StsGroupBox();
    JTextField jTextField1 = new JTextField();

    public StsPowerPanel()
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
        jLabel1.setText("Power");
        this.setLayout(gridLayout1);
        groupBox1.setLabel("Inverse Distance Parameters");
        jTextField1.setText("2                 ");
        this.add(groupBox1, null);
        groupBox1.add(jLabel1, null);
        groupBox1.add(jTextField1, null);
    }

    public float getPower() { return Float.valueOf(jTextField1.getText()).floatValue(); }

    public static void main(String[] args)
    {
        StsPowerPanel powerPanel = new StsPowerPanel();
		com.Sts.Utilities.StsToolkit.createDialog(powerPanel);
    }
}
