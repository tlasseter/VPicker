package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Wizards.PlatformPlan.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
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

public class StsSelectPlatformSlot extends StsWizardStep
{
    StsConfigurePlatformPanel panel;
    StsHeaderPanel header;

    StsPlatform platform = null;

    StsDefineWellPlanWizard wizard = null;

    public StsSelectPlatformSlot(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsDefineWellPlanWizard)wizard;
        panel = new StsConfigurePlatformPanel(wizard);
        panel.setEditable(false);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 700));
        header.setTitle("Select A Platform Slot for this Planned Well");
        header.setSubtitle("Interactively pick slot on diagram");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#WellPlan");
        header.setInfoText(wizardDialog,"Left click on a slot which doesn't already have a well assigned.");
    }

    public boolean start()
    {
        platform = (StsPlatform)wizard.getPlatform();
        panel.initialize(platform);
        wizard.rebuild();
        return true;
    }
/*
    public boolean cancel()
    {
        platform.deleteWellAtCurrentSlot();
        return false;
    }
*/
    public boolean end()
    {
        if(!wizard.isCanceled())
        {
            if (platform.getCurrentSlotIndex() == -1) {
                new StsMessage(wizard.getModel().win3d, StsMessage.WARNING,
                               "No slot has been selected.");
                return false;
            }
            double[] slotXY = platform.getSlotXY();
            if(!wizard.addWellToPlatform()) return true;
            wizard.setWellXY(slotXY);
            wizard.setZKB(platform.getZKB());
            return true;
        }
        return false;
    }
}
