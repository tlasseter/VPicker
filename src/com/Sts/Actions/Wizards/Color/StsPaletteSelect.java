
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Color;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.IO.*;

import java.awt.*;

public class StsPaletteSelect extends StsWizardStep
{
    StsPaletteSelectPanel panel;
    StsHeaderPanel header;

    public StsPaletteSelect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsPaletteSelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 300));
        header.setTitle("Palette Selection");
        header.setSubtitle("Selecting Available Palettes");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Color");
        header.setInfoText(wizardDialog,"(1) Navigate to the directory containing the Palette Files using the" +
                           " Dir Button. All palette files in the selected directory will be placed" +
                           " in the left list.\n      ***** File names must have a palette.txt. prefix *****\n" +
                           "(2) Select the desired files from the left list and place them in the" +
                           " right list using the provided controls between the lists.\n(3) Palettes can be" +
                           " viewed on the left side by selecting the desired palette from either list\n" +
                           "(4) Once all desired palettes have been selected press the Next>> Button");
    }

    public StsAbstractFile[] getSelectedFiles() {  return panel.getSelectedFiles(); }
    public String[] getSelectedPaletteNames() {  return panel.getSelectedPaletteNames(); }
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

