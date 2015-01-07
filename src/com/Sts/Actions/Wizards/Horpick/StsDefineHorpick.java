package com.Sts.Actions.Wizards.Horpick;

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

public class StsDefineHorpick extends StsWizardStep
{
    StsDefineHorpickPanel panel;
    StsHeaderPanel header;
    StsHorpick horpick;

    public StsDefineHorpick(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineHorpickPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 600));
        header.setTitle("Pick Horizon Definition");
        header.setSubtitle("Defining Horizon for Picking");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#HorPick");        
        header.setInfoText(wizardDialog,"(1) Specify the desired horizon name.\n" +
                                  "(2) Select the surface color (may change from object panel after creation).\n" +
                                  "(3) Enter the Pick Preferences Dialog to adjust graphical refresh rate during picking.\n" +
                                  "   ***** The defaults rate is every pick cycle and should not be adjusted unless ***** \n" +
                                  "   ***** graphics refresh appears to be delayed. *****\n" +
                                  "(4) Press the Next>> Button to proceed to interactive picking.");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/HorizonTracking.html");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
//        return constructHorpick();
        return true;
    }

    private boolean constructHorpick()
    {
        try
        {
            String name = panel.getName();
            if(name.length() == 0)
            {
                new StsMessage(wizard.frame, StsMessage.ERROR, "Please enter a name for this horizon.");
                return false;
            }

            horpick = (StsHorpick)model.getObjectWithName(StsHorpick.class, name);
            if(horpick == null)
                horpick = StsHorpick.constructor(name, panel.getStsColor(), (StsHorpickWizard)wizard);
            if(horpick == null) return false;

            model.setCurrentObject(horpick);
            return true;
        }
        catch(Exception e)
        {
            new StsMessage(wizard.frame, StsMessage.ERROR, "Failed to construct Horpick.\n" + e.getMessage());
            return false;
        }
    }
    
    public StsHorpick getDefinedHorpick()
    {
        if(!constructHorpick()) return null;
        return horpick;
    }
}
