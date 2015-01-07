package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;

public class StsSensorDefineStages extends StsWizardStep {
    StsSensorDefineStagesPanel panel;
    StsHeaderPanel header;

    public StsSensorDefineStages(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSensorDefineStagesPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Sensor Load");
        header.setSubtitle("Define Sensor Colors");
        header.setInfoText(wizardDialog, "(1) Specify the Number of Stages\n" +
                                        "(2) Verify and Adjust Stage Times\n" +
                                        "(3) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SensorLoad");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }

}