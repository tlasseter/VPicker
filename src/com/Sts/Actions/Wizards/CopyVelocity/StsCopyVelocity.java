package com.Sts.Actions.Wizards.CopyVelocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsCopyVelocity extends StsWizardStep
{
    StsCopyVelocityPanel panel;
    StsHeaderPanel header;

    public StsCopyVelocity(StsWizard wizard)
    {
        super(wizard);
        panel = new StsCopyVelocityPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 500));
        header.setTitle("Copy Velocity Model");
        header.setSubtitle("Define Velocity Model");
        header.setInfoText(wizard.dialog,"(1) Specify the new velocity model name.\n" +
                                         "(2) Press the Finish Button to create a copy.");
    }

    public boolean start()
    {
        wizard.enableFinish();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }
}
