
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsPreStackSelect extends StsPostStackSelect
{
    public StsPreStackSelect(StsPreStackWizard wizard)
    {
        super(wizard);
        panel = new StsPostStackSelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Pre-Stack SegY Definition");
        header.setSubtitle("Select Pre-Stack3d Files(s)");
//        panel.setPreferredSize(new Dimension(500, 600));
        header.setInfoText(wizardDialog,"(1) Navigate to the directory containing the SegY Files using the Dir Button.\n" +
                                        "    All SegY Files in the selected directory will be placed in the left list.\n" +
                                        "      **** File names must have a .sgy, .SGY, .segy or .SEGY suffix *****\n" +
                                        "(2) Select files from the left list and place them in the right using the controls between the lists.\n" +
                                        "      **** All selected files will be scanned and placed in the table at bottom of screen with stats****\n" +
                                        "      **** from a random scan for review to ensure files are read correctly and disk space is adequate. ****\n" +
                                        "      **** If scan fails, adjust parameters on this and next screen before proceeding w/ processing ****\n" +
                                        "(3) Select a pre-defined Segy Format template if one exists for the selected files\n" +
                                        "      **** Format templates are saved from previous executions of this wizard.\n" +
                                        "(4) Set the scan percentage (# of samples of the SegY Files scanned) for verifying correct file reading.\n" +
                                        "      **** Random scanning of samples is used to determine if file is being read correctly.\n" +
                                        "(5) Once all selected file formats are correctly defined press the Next>> Button\n" +
                                        "      **** Dont worry if value in table are incorrect, next screen allows trace header mapping.****");
    }
}

