
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.VolumeStimulated;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsSelectFractureSets extends StsWizardStep
{
    StsSelectFractureSetsPanel panel;
    StsHeaderPanel header;

    public StsSelectFractureSets(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectFractureSetsPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Select Fracture Sets");
        header.setSubtitle("Select the Fracture Sets for Analysis");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VolumeStimulated");
        header.setInfoText(wizardDialog,"(1) Select a fracture set from the list.\n" +
        		"(2) Select the Enable switch to include or exclude selected set from analysis.\n" +
        		"(3) Specify the radius and cell size of the stimulated area for the selected set.\n" +
        		"(4) Repeat for all desired fracture sets, confirm coverage in 3D view.\n" +
        		"(5) Press the Run button to compute the volume stimulated.\n" +
        		"(6) Selected units for volume calculation and re-Run to compute in selected units.\n" +
        		"(7) Press the Next>> Button");
    }

    public boolean start()
    {
    	enableFinish();
        wizard.disableNext();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }
    
    public void setAnalysisMessage(String msg)
    {
    	panel.setMessage(msg);
    }    
}

