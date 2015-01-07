package com.Sts.Actions.Wizards.SimulationFile;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.RFUI.*;

import java.awt.*;
import java.io.*;

public class StsSelectSimulationFile extends StsWizardStep
{
    StsSimulationFileChoicePanel panel;
    StsHeaderPanel header;

    public StsSelectSimulationFile(StsWizard wizard)
    {
        super(wizard, new StsSimulationFileChoicePanel(), null, new StsHeaderPanel());
        panel = (StsSimulationFileChoicePanel) getContainer();
        setDefaultDirectory(model.getProject().getRootDirectory());
        panel.setPreferredSize(new Dimension(400,300));

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Simulation File Export");
        header.setSubtitle("Select Simulation File");
    }

    public boolean start()
    {
        return true;
    }


    public boolean end() { return getFilename() != null; }
    public int getFormat() { return panel.getFormat(); }
    public String getFilename() { return panel.getFilename(); }
    public void setDefaultDirectory(File dir) { panel.setCurrentDir(dir); }

}
