
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Horpick;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.StsSeismicVolume;
import com.Sts.DBTypes.StsSurface;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsYesNoDialog;

import java.awt.*;

public class StsSelectHorpick extends StsWizardStep
{
    StsSelectHorpickPanel panel;
    StsHeaderPanel header;

    public StsSelectHorpick(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectHorpickPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 600));
        header.setTitle("Horizon Selection");
        header.setSubtitle("Selecting Available Horizons");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#HorPick");                
        header.setInfoText(wizardDialog,"(1) Select the surface that you wish to edit or press the New Surface button.\n" +
                                  "   ***** If editing an existing surface, all previously picked seeds ***** \n" +
                                  "   ***** and picks are retained and can be refined, and re-picked *****\n" +
                                  "(2) Press the Next>> Button if selecting an existing surface.");

        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/HorizonTracking.html");
    }

    public boolean start()
    {
        // Can only pick surfaces in the native domain of the seismic data so need to confirm
        // volumes exist in current domain and if not switch domains.
        StsSeismicVolume[] volumes = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
        StsSeismicVolume vol = volumes[0];
        volumes = StsSeismicVolume.getZDomainVolumes(volumes, model.getProject().getZDomain());
        if(volumes.length == 0)
        {
            if(StsYesNoDialog.questionValue(wizard.frame, "No volumes exist in " + model.getProject().getZDomainString() + ".\nDo you want to switch domains?"))
                model.getProject().setZDomain(vol.getZDomain());
            else
            {
                new StsMessage(wizard.frame, StsMessage.INFO, "Unable to pick surface in non-native domain. Exiting wizard.");
                wizard.cancel();
                return false;
            }
        }
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }
}

