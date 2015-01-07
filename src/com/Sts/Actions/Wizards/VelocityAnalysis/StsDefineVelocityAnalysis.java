
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.VelocityAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsDefineVelocityAnalysis extends StsWizardStep
{
    StsDefineVelocityAnalysisPanel panel;
    StsHeaderPanel header;

    public StsDefineVelocityAnalysis(StsWizard wizard)
    {
        super(wizard);
        boolean is3d = StsPreStackLineSetClass.currentProjectPreStackLineSet instanceof StsPreStackLineSet3d;
        if (is3d) {
        	panel = new StsDefineVelocityAnalysisPanel3d(wizard, this);
        } else {
        	panel = new StsDefineVelocityAnalysisPanel2d(wizard, this);
        }
        	
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 500));
        header.setTitle("Velocity Analysis Definition");
        header.setSubtitle("Analysis Point Setup");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VelocityAnalysis");
        header.setInfoText(wizard.dialog,"(1) Specify the row and column (CDP) increment for velocity analysis.\n" +
                           "      **** Increments will be used to set-up an analysis grid used to define\n" +
                           "      **** and maintain status of velocity analysis ****\n" +
                           "(2) Specify the trace threshold required for velocity analysis\n." +
                           "      **** Trace threshold is the minimum number of traces required to allow a \n" +
                           "      **** velocity profile to be picked. ****\n" +
                           "(3) Press the Next>> Button.");
    }

    public boolean start()
    {
        //wizard.enableFinish();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }
}

