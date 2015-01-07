
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Pointset;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.IO.*;

import java.awt.*;

public class StsPointsetSelect extends StsWizardStep
{
    StsPointsetSelectPanel panel;
    StsHeaderPanel header;

    public StsPointsetSelect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsPointsetSelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 300));
        header.setTitle("Pointset Selection");
        header.setSubtitle("Selecting Available Pointsets");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/PointsetLoad.html");
        header.setInfoText(wizardDialog,"(1) Navigate to the directory containing the Pointset Files using the" +
                           " Dir Button. All pointset files in the selected directory will be placed" +
                           " in the left list.\n      ***** File names must have a .csv prefix *****\n" +
                           "(2) Select the desired files from the left list and place them in the" +
                           " right list using the provided controls between the lists.\n" +
                           "(3) Once all desired palettes have been selected press the Next>> Button");
    }

    public StsAbstractFile[] getSelectedFiles() {  return panel.getSelectedFiles(); }
    public String[] getSelectedPointsetNames() {  return panel.getSelectedPointsetNames(); }
    public String getSelectedDirectory() {  return panel.getSelectedDirectory(); }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }
}

