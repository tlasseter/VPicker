package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsColumnDefine extends StsWizardStep
{
    StsColumnDefinePanel panel;
    StsHeaderPanel header;

    public StsColumnDefine(StsWizard wizard)
    {
        super(wizard);
        panel = new StsColumnDefinePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
//        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Well Load");
        header.setSubtitle("Column Definition");
        header.setInfoText(wizardDialog,"(1) \n" +
                           "(2) Press the Next>> Button.\n");
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

}
