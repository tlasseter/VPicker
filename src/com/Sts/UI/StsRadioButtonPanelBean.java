
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;


import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;

class TestPanelObject
{
    public byte buttonState;
    public String[] buttonLabels = new String[] { "One", "Two", "Three" };

    public TestPanelObject()
    {
    }

    public String getButtonState() { return buttonLabels[buttonState]; }
    public void setButtonState(String buttonStateString)
    {
        for(int n = 0; n < 3; n++)
            if(buttonStateString == buttonLabels[n])
            {
                buttonState = (byte)n;
                System.out.println("Radio button " + buttonStateString + " selected.");
                return;
            }
        System.out.println("Failed to find button selected.");
    }
}

