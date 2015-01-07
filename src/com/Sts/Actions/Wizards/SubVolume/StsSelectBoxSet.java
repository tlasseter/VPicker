
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SubVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsSelectBoxSet extends StsWizardStep
{
    StsSelectBoxSetPanel panel;
    StsHeaderPanel header;
    StsBoxSetSubVolume selectedBoxSet = null;

    public StsSelectBoxSet(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectBoxSetPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Sub-Volume Definition");
        header.setSubtitle("Selecting Box Set");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SubVolume");
        header.setInfoText(wizardDialog, "(1) Select an existing Box Set Subvolume or press the New SubVolume Button.\n" +
                           "     **** Only Box Set SubVolumes are editable from the wizard ****\n" +
                           "     **** Other SubVolume types can be edited from the object panel ****\n" +
                           "(2) Press the Next> Button once the desired Box Set is selected.\n");
    }

    public boolean start()
    {
        panel.initialize();
		StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
		subVolumeClass.setIsApplied(false);
        return true;
    }

    public boolean end()
    {
//		StsBoxSetSubVolumeClass boxSetClass = getBoxSetClass();
//		boxSetClass.setIsApplied(false);
//		boxSetClass.setIsVisible(true);
        return true;
    }

    public void setBoxSet(StsBoxSetSubVolume boxSet)
    {
		if(selectedBoxSet != null) selectedBoxSet.setIsVisible(false);
		this.selectedBoxSet = boxSet;
		selectedBoxSet.setIsVisible(true);
        model.setCurrentObject(boxSet);
    }

    public StsBoxSetSubVolume getBoxSet()
    {
        return selectedBoxSet;
    }

    // if boxSet is null, defineBoxSetSubVolume will create a new one
    public void createBoxSet()
    {
//        selectedBoxSet = null;
        model.getStsClass(StsBoxSetSubVolume.class).setCurrentObject(null);
        wizard.next();
    }
}

