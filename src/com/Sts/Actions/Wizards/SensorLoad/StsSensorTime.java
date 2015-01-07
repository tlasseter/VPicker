package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsSensorTime extends StsWizardStep {
    StsSensorTimePanel panel;
    StsHeaderPanel header;

    public StsSensorTime(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSensorTimePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Sensor Load");
        header.setSubtitle("Define Sensor Time");
        header.setInfoText(wizardDialog, "(1) Enter the starting time.\n" +
        		           " **** Make sure the month-day order agrees with the Date Order in the project.\n" +
                           " **** If supplied time is relative, only the supplied date will be used.\n" +
                           "(2) Press the Next>> Button.\n");
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
