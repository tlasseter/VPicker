
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SimulationFile;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.Horizons.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Actions.Export.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;

public class StsDefineMinCellSize extends StsWizardStep
{
    StsJPanel panel;
    StsHeaderPanel header;

    public StsDefineMinCellSize(StsWizard wizard)
    {
    	super(wizard);
        panel = StsEclipseOutput.getEclipseMergeCellPanel(wizard.model);
        this.setPanels(panel, new StsHeaderPanel());
        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Cell Truncation Specification");
        // header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Horizons");
//        panel.setPreferredSize(new Dimension(400,300));
        header.setInfoText(wizardDialog,"(1) Select the minimum allowable cell truncation fraction.\n" +
                                        "    The number of merged cells and the number of cells smaller than this fraction will be computed.\n " +
                                        "(2) When satisfied, Press the Next>> Button.");
    }

    public boolean start()
    {
        return true;
    }

    public void cancel()
    {
        StsEclipseOutput.cancel();
    }

    public boolean end()
    {
        return true;
    }
}