package com.Sts.Actions.Wizards.VolumeStimulated;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.Well.StsFileTypePanel;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;

public class StsSelectType extends StsWizardStep {
    StsSelectTypePanel panel;
    StsHeaderPanel header;

    public StsSelectType(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectTypePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Volume Stimulated");
        header.setSubtitle("Select Computation Type");
        header.setInfoText(wizardDialog,"(1) Select how to compute the stimulated volume.\n" +
                            "       Fracture set will compute from a radius defined around selected\n" +
                            "       fracture sets.\n" +
                            "       Gridded sensor volume will compute from filled cells in the selected\n" +
                            "       sensor volume.\n" +
                           "(2) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VolumeStimulated");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public byte getType()
    {
        return panel.getType();
    }

}