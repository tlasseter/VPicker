package com.Sts.Actions.Wizards.Culture;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsCultureSetAttributes extends StsWizardStep
{
	StsCultureSetAttributesPanel panel;
    StsHeaderPanel header;
    StsHorpick horpick;

    public StsCultureSetAttributes(StsWizard wizard)
    {
        super(wizard);
        panel = new StsCultureSetAttributesPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 600));
        header.setTitle("Culture Set Definition");
        header.setSubtitle("Defining Imported Cultural Set");
        header.setInfoText(wizardDialog,"(1)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Culture");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }

}
