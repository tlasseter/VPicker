
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

public class StsEditFaultGaps extends StsWizardStep
{
    StsEditFaultGapsPanel panel;
    StsHeaderPanel header;

    public StsEditFaultGaps(StsWizard wizard)
    {
        super(wizard, new StsEditFaultGapsPanel(), null, new StsHeaderPanel());
        panel = (StsEditFaultGapsPanel) getContainer();
        panel.setTitle("Modify fault gaps before building the 3d grid:");
        panel.setPreferredSize(new Dimension(400,300));

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Model Construction");
        header.setSubtitle("Edit Fault Gaps");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Model");                
        header.setInfoText(wizardDialog,"(1) Select a fault section from the list on left. \n" +
                                 "(2) Specify the number of cells to remove around fault during model construction.\n" +
                                 "   ***** If horizons were picked in S2S, gap back can be made *****\n" +
                                 "   ***** by pick quality on the next panel.  ***** \n\n" +
                                 "   ***** Gaps can be defined for all fault sections at once using multi-select. *****\n" +
                                 "(3) Press the Next>> Button");
    }

    public boolean start()
    {
    	panel.setSections(model.getObjectListOfType(StsSection.class, StsSection.FAULT));
        return true;
    }

    public boolean end()
    {
        return true;
    }
}
