package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsSensorRelative extends StsWizardStep {
    StsSensorRelativePanel panel;
    StsHeaderPanel header;

    public StsSensorRelative(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSensorRelativePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Sensor Load");
        header.setSubtitle("Define Absolute Location");
        header.setInfoText(wizardDialog, "**** Some of the selected files contain relative coordinates ****\n" +
                            "(1) Specify the absolute coordinates to add to the relative coordinates in the file\n" +
                            " **** 0.0, 0.0, 0.0 position will leave the sensors in the positions supplied in the file ****\n" +
                           "(2) Press the Next>> Button.");
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
