package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsPlatformOrWell extends StsWizardStep
{
    StsPlatformOrWellPanel panel;
    StsHeaderPanel header;

    StsDefineWellPlanWizard wizard = null;

    StsWellPlan plan;
    StsWellPlan planSet = null;
    StsCursor3d cursor3d;
    StsSeismicVelocityModel velocityModel;

    StsColor stsColor;
    String name = new String("WellPlan");
    Color color = Color.RED;
    boolean visible = false;

    static int suffix = 1;

    public StsPlatformOrWell(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsDefineWellPlanWizard)wizard;
        panel = new StsPlatformOrWellPanel(wizard);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Define name and well Location Type");
        header.setSubtitle("Platform or pad, pattern, or single well");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#WellPlan");
        header.setInfoText(wizardDialog, "Set name and type and then go to next step.\n");
    }

    public boolean start()
    {
//        planSet = wizard.getCreatePlanSet();
//        panel.classInitialize(planSet);
		panel.initialize();
        return true;
    }

    public void cancel()
    {
        wizard.setWellPlanSet(null);
    }

    public boolean end()
    {
        return true;
    }
}
