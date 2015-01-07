package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;

public class StsWellHeadDefinition extends StsWizardStep {
    StsWellHeadDefinitionPanel panel;
    StsHeaderPanel header;

    public StsWellHeadDefinition(StsWizard wizard)
    {
        super(wizard);
        panel = new StsWellHeadDefinitionPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Well Load");
        header.setSubtitle("Wellhead Definition");
        header.setInfoText(wizardDialog,"(1) Select the desired well(s)\n" +
                            "(2) Specify the coordinates of the wellhead. Defaults are read from selected file.\n" +
                            "(3) Select whether to subtract the kelly bushing elevation from supplied depths\n" +
                            "(4) Override the kelly bushing elevation supplied in the selected files\n" +
                           "(5) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Well");
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

    public int[] getSelectedIndices()
    {
       return panel.getSelectedIndices();
    }

}