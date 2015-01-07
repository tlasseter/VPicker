
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PreStack3dLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.IO.*;

public class StsPreStackLoadSelect extends StsWizardStep
{
    StsPreStackLoadSelectPanel panel;
    StsHeaderPanel header;

    public StsPreStackLoadSelect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsPreStackLoadSelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
//        panel.setPreferredSize(new Dimension(400, 450));
        header.setTitle("Pre-Stack Line Selection");
        header.setSubtitle("Select Available Lines");
        header.setInfoText(wizardDialog,"(1) Select the seismic line(s) that will be loaded.\n" +
                           "     ***** Must pre-process SegY Files before loading them. See Process Pre-Stack Seismic step. *****\n" +
                           "(2) If other lines are desired, navigate to the directory containing the lines using the Dir button.\n" +
                           "     ***** All lines in the selected directory will be placed in the left list. *****\n" +
                           "     ***** The default location is the project directory but any directory can be selected ****\n" +
                           "(3) Select the desired lines from the left list and place them in the right list using the provided controls between the lists.\n" +
                           "     ***** Information regarding the currently selected line is isVisible in the panel at the bottom of the screen.*****\n" +
                           "(4) An initial palette can be selected using the Select Palette Button or left to the default Red-White-Blue palette.\n" +
                           "(5) Once all line selections are complete, press the Next>> Button.");
    }

    public StsFile[] getSelectedFiles()
    {
        return panel.getSelectedFiles();
    }

    public String[] getFilenameEndings(StsAbstractFile[] files)
    {
        return panel.getFilenameEndings(files);
    }

    public boolean start()
    {
        wizard.dialog.setTitle("Select Available Lines");
        return true;
    }

    public boolean end()
    {
        return true;
    }
}

