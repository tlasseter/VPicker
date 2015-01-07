
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SurfacesFromMarkers;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsSelectMarkers extends StsWizardStep
{
    StsSelectMarkersPanel panel;
    StsHeaderPanel header;

    public StsSelectMarkers(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectMarkersPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 300));
        header.setTitle("Create Surfaces from Markers");
        header.setSubtitle("Select Markers");
        header.setInfoText(wizard.dialog, "(1) Select one or more markers from the list.\n" +
            "    **** Multi-select using the shift or control and left mouse button combination\n" +
            "(2) Specify the number of inline(Y) and crossline(X) grid cells or bins\n" +
            "(3) Specify the grid origin\n" +
            "(4) Specify the X and Y interval in feet or meters."
            + "(5) Press the Next Button.");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SurfacesFromMarkers");
    }

    public StsMarker[] getSelectedMarkers()
    {
        return panel.getSelectedMarkers();
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }
}

