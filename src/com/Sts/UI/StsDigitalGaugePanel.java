
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;


import com.Sts.UI.Beans.StsJPanel;
import com.nextwavesoft.enumeration.*;
import com.nextwavesoft.gauge.*;
import com.nextwavesoft.primitives.GaugeItem;
import com.nextwavesoft.shared.Unit;

import javax.swing.*;
import java.awt.*;

public class StsDigitalGaugePanel extends StsJPanel
{
    StsDigitalGauge digitalGauge = null;

	public StsDigitalGaugePanel()
    {
		try
        {
            createGauge();
        }
        catch(Exception e) { e.printStackTrace(); }
    }

    private void createGauge()
    {

        gbc.fill = gbc.BOTH;
        digitalGauge = new StsDigitalGauge();
        addEndRow(digitalGauge);
    }

    public void setValue(float value)
    {
        digitalGauge.setText(Float.toString(value));
    }

    public static void main(String[] argv)
    {
        float min = 0.0f;
        float max = 1000.0f;

        //Create and set up the window.
		JFrame frame = new JFrame("Digital Gauge");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        StsDigitalGaugePanel m = new StsDigitalGaugePanel();
        frame.setContentPane(m);

		//Display the window.
		frame.pack();
		frame.setSize(new Dimension(200, 80));
		frame.setVisible(true);

        float inc = (max-min)/1000;
        for(float i=min; i<=max; i += inc)
        {
            try { Thread.sleep(10); }
            catch(Exception ex) {}
            m.setValue(i);
        }
    }
}