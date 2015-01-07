package com.Sts.Actions.Wizards.Horpick;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsFinishHorpick extends StsWizardStep
{
    public StsFinishHorpickPanel panel;
    public StsHeaderPanel header;

    public StsFinishHorpick(StsWizard wizard)
    {
        super(wizard);
        panel = new StsFinishHorpickPanel(wizard);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 600));
        header.setTitle("Complete Horizon Definition");
        header.setSubtitle("Smooth and Fill");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#HorPick");                
        header.setInfoText(wizardDialog,"(1) Quality control the surface by applying attributes to the surface\n" +
                                  "   **** Look for cycle skips in colors. ****\n" +
                                  "(2) Apply the correlation coefficeints to the surface.\n" +
                                  "   **** Adjust color spectrum as desired ****\n" +
                                  "(3) Using the horizontal slider, remove poor quality picks and mis-picks.\n" +
                                  "   **** If picks are disconnected from the seed that generate them during **** \n" +
                                  "   **** this process, all disconnectted picks will also be removed. ****\n" +
                                  "(4) Press the Next>> Button to finish the surface picking\n" +
                                  "   **** The surface picker can be re-entered and picks, seeds ****\n" +
                                  "   **** and surface adjusted at any time.");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/HorizonTracking.html");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
//        return wizard.end();
    }
}
