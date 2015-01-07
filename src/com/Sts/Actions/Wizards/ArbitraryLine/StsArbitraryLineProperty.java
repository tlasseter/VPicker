package com.Sts.Actions.Wizards.ArbitraryLine;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsArbitraryLineProperty extends StsWizardStep
{
    StsArbitraryLinePropertyPanel panel;
    StsHeaderPanel header;

    public StsArbitraryLineProperty(StsWizard wizard)
    {
        super(wizard);
        panel = new StsArbitraryLinePropertyPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 400));
        header.setTitle("Arbitrary Line Definition");
        header.setSubtitle("Defining Line Properties");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ArbitraryLine");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        try
        {
            String name = panel.getName();
            if(name.length() == 0)
            {
                new StsMessage(wizard.frame, StsMessage.ERROR, "Please enter a name for this arbitrary line.");
                return false;
            }
        }
        catch(Exception e)
        {
            new StsMessage(wizard.frame, StsMessage.ERROR, "Failed to construct line.\n" + e.getMessage());
            return false;
        }
        return true;
    }
}
