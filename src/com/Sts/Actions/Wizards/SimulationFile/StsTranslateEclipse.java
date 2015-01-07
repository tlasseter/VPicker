package com.Sts.Actions.Wizards.SimulationFile;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.RFUI.*;

public class StsTranslateEclipse extends StsWizardStep
{
    StsKxyzPanel panel;
    StsHeaderPanel header;

    public StsTranslateEclipse(StsWizard wizard)
    {
        super(wizard);
        dialogPanel = new StsTranslateEclipsePanel(model);
        headerPanel = new StsHeaderPanel(); 
        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Translate Eclipse <-> S2S Indexes & Coordinates");
        header.setSubtitle("");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }
}