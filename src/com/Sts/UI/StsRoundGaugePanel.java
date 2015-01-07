
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;


import com.Sts.DBTypes.StsObject;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.Utilities.StsParameters;
import com.nextwavesoft.enumeration.*;
import com.nextwavesoft.gauge.*;
import com.nextwavesoft.primitives.GaugeItem;
import com.nextwavesoft.shared.Unit;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ItemListener;
import java.util.ArrayList;

public class StsRoundGaugePanel extends StsJPanel
{
    StsRoundGauge gauge = null;

    public StsRoundGaugePanel()
    {
        this(null, 0, 200, 0, new Color(44,44,44));
    }

	public StsRoundGaugePanel(String title, float min, float max)
    {
        this(title, min, max, 0f, new Color(44,44,44));
    }

    public StsRoundGaugePanel(String title, float min, float max, float highlightMin, Color rim)
    {
		try
        {
            Color rimColor = new Color(44,44,44);
            if(rim != null)
                rimColor = rim;
            createGauge(title, min, max, highlightMin, rimColor);
        }
        catch(Exception e) { e.printStackTrace(); }
    }


    public StsRoundGaugePanel(String title, float min, float max, int highlightPercentage, Color rim)
    {
		try
        {
            Color rimColor = new Color(44,44,44);
            float highlightMin = 0f;
            if(highlightPercentage != 0)
                highlightMin = max - ((max-min)*(1.f - (highlightPercentage/100f)));
            if(rim != null)
                rimColor = rim;
            createGauge(title, min, max, highlightMin, rimColor);
        }
        catch(Exception e) { e.printStackTrace(); }
    }

    private void createGauge(String title, float min, float max, float highlightMin, Color rimColor)
    {
        gauge = new StsRoundGauge(title, min, max, highlightMin, rimColor);
        gbc.fill = gbc.BOTH;
        addEndRow(gauge);
    }

    public void setValue(float value)
    {
        if(value == StsParameters.nullValue)
        {
            gauge.setValue(gauge.getMin());
            gauge.setDigitalValue("----.-");
        }
        else
            gauge.setValue(value);
    }

    public static void main(String[] argv)
    {
        float min = 0f;
        float max = 1000f;

        //Create and set up the window.
		JFrame frame = new JFrame("RoundGauge Display");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        StsRoundGaugePanel m = new StsRoundGaugePanel("magnitude", min, max, 90, new Color(0,0,180));
        frame.setContentPane(m);

		//Display the window.
		frame.pack();
		frame.setSize(new Dimension(400, 400));
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