
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Model;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsSurfaceAutopick extends StsWizardStep
{
    StsSurfaceAutopickPanel panel;
    StsHeaderPanel header;

    public StsSurfaceAutopick(StsWizard wizard)
    {
		super(wizard, new StsSurfaceAutopickPanel(wizard.getModel()), null, new StsHeaderPanel());
        panel = (StsSurfaceAutopickPanel)dialogPanel;
        panel.setTitle("Modify units before building the 3d grid:");
        panel.setPreferredSize(new Dimension(400,300));

		header = (StsHeaderPanel)headerPanel;
        header.setTitle("Model Construction");
        header.setSubtitle("Specify Fault Gaps base on Pick Quality");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Model");                
        header.setInfoText(wizardDialog,"(1) Select a horizon from the left list.\n" +
                           "   ***** Horizon must have been picked in S2S, otherwise skip this screen *****\n" +
                                 "(2) Specify the pick quality to accept in model construction.\n" +
                                 "   ***** All picks below specified value will be removed and interpolated through. *****\n" +
                                 "(3) Press the Next>> Button");
    }

    public boolean start()
    {
        StsModelSurface[] surfaces =(StsModelSurface[])model.getCastObjectList(StsModelSurface.class);
    	panel.setSurfaces(surfaces);
        return true;
    }

    public boolean end()
    {
        return true;
    }
}
