
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Model;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsEditZones extends StsWizardStep
{
    StsEditZonesPanel panel;
    StsHeaderPanel header;

    public StsEditZones(StsWizard wizard)
    {
        super(wizard, new StsEditZonesPanel(), null, new StsHeaderPanel());
        panel = (StsEditZonesPanel) getContainer();
        panel.setTitle("Modify units before building the 3d grid:");
        panel.setPreferredSize(new Dimension(400,300));

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Model Construction");
        header.setSubtitle("Define Zones");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Model");                
        header.setInfoText(wizardDialog,"(1) Select a zone from the list on left.\n" +
                                 "(2) Change the name, type and number of sub-zones.\n" +
                                 "   ***** Number of sub-zones is the initial sub-zone numbers *****\n" +
                                 "   ***** and can be modified after model construction. ***** \n" +
                                 "(3) Press the Next>> Button");
    }

    public boolean start()
    {
        StsZone[] zones =(StsZone[])model.getCastObjectList(StsZone.class);
    	panel.setZones(zones);
        return true;
    }

    public boolean end()
    {
        return true;
    }
}
