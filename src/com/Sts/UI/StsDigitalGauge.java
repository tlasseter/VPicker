
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;


import com.Sts.UI.Beans.StsJPanel;
import com.nextwavesoft.enumeration.DigitalCharacterOverflowAlignment;
import com.nextwavesoft.enumeration.DigitalCharacterType;
import com.nextwavesoft.enumeration.RectangularFrameType;
import com.nextwavesoft.gauge.DigitalGauge;

import javax.swing.*;
import java.awt.*;

public class StsDigitalGauge extends DigitalGauge
{
	public StsDigitalGauge()
    {
        super();
        configureDigitalMeter();
    }

    /**
     * Configure a digital meter.
     *
     * @return the digital gauge
     */
    public void configureDigitalMeter()
    {
		setFrameType(RectangularFrameType.RoundedSideThinRim);
		setBackground(Color.black);
		setCharacterType(DigitalCharacterType.Segment7Trapezoid);
        setCharacterCount(10);
		setForeground(Color.WHITE);
		setDimmedBrush(Color.black);
		setCharacterOverflowAlignment(DigitalCharacterOverflowAlignment.Right);
        setCharacterPadding(1);
        setCharacterSpacing(1);
		setItalic(false);
    }

}