package com.Sts.Actions.Wizards.PlatformPlan;

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

public class StsConfigurePlatform extends StsWizardStep
{
    StsConfigurePlatformPanel panel;
    StsHeaderPanel header;

    StsPlatform platform = null;

    StsPlatformPlanWizard wizard = null;

    public StsConfigurePlatform(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsPlatformPlanWizard)wizard;
        panel = new StsConfigurePlatformPanel(wizard);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 700));
        header.setTitle("Platform Plan");
        header.setSubtitle("Configure Platform");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#PlatformPlan");
        header.setInfoText(wizardDialog,"(1) Specify the name of the plan.\n");
    }

    public boolean start()
    {
        platform = wizard.getPlatform();
        panel.initialize(platform);
        return true;
    }

    public boolean end()
    {
        platform.configurePlatform();
        return true;
    }

}
