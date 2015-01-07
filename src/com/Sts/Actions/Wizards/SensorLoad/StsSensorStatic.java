package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsSensorStatic extends StsWizardStep {
    StsSensorStaticPanel panel;
    StsHeaderPanel header;

    public StsSensorStatic(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSensorStaticPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Sensor Load");
        header.setSubtitle("Define Static Location");
        header.setInfoText(wizardDialog, "**** Some files do not contain X, Y, Z values so static position is required ****\n" +
                            "(1) Specify the absolute coordinates of the sensor. \n" +
                            " **** Static sensors are initialized s invisible in 2D-3D view ****\n" +
                            " **** they can be enabled from the Object Panel if desired. ****\n" +
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
