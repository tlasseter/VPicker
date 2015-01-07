package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.*;

public class StsHeaderFile extends StsWizardStep {
    StsHeaderFilePanel panel;
    StsHeaderPanel header;

    public StsHeaderFile(StsWizard wizard)
    {
        super(wizard);
        panel = new StsHeaderFilePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
//        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Well Load");
        header.setSubtitle("Header File Definition");
        header.setInfoText(wizardDialog,"(1) \n" +
                           "(2) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Well");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        if(panel.getHeaderFile() == null)
        {
            new StsMessage(wizard.frame, StsMessage.WARNING, "Header file selection is required to proceed.");
            return false;
        }
        return true;
    }

}
