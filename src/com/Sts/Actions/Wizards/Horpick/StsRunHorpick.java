package com.Sts.Actions.Wizards.Horpick;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsRunHorpick extends StsWizardStep
{
    public StsRunHorpickPanel panel;
    public StsHeaderPanel header;

    public StsRunHorpick(StsWizard wizard)
    {
        super(wizard);
        panel = new StsRunHorpickPanel(wizard);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 600));
        header.setTitle("Pick Horizon Definition");
        header.setSubtitle("Defining Horizon for Picking");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#HorPick");                
        header.setInfoText(wizardDialog,"(1) Select the mode of operation.\n" +
                                  "   **** Seed - Select seeds in the graphics window. ***** \n" +
                                  "   **** Delete Patch - Delete the current patch and all associated picks ****\n" +
                                  "   **** Delete Picks - Delete the picks of the current patch, retain the seed point. ****\n" +
                                  "   **** The current patch is the one selected in the dropdown list ****\n" +
                                  "(2) Specify the event, window length and pick difference for the seed point prior to picking it.\n" +
                                  "   **** Event is maximum, minimum, positive zero cross or negative zero cross. ****\n" +
                                  "   **** Window Length is the correlation window to be tracked. ****\n" +
                                  "   **** Pick Difference is effectively the allowable dip. ****\n" +
                                  "(3) Select the seed point and review the seed properties.\n" +
                                  "   **** Make numerous seed picks anywhere discontinuities in the event occur. ****\n" +
                                  "(4) Track the seed picks individually or all by pressing the Run button.\n" +
                                  "   **** Interative will run all seeds with highest criteria and then step down and run again. ****\n" +
                                  "(5) Press the Next>> Button when surface is propoerly tracked.");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/HorizonTracking.html");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public void setStatusMessage(String text)
    {
        panel.setStatusLabel(text);
    }

    public boolean end()
    {
//        return wizard.end();
        return true;
    }
}
