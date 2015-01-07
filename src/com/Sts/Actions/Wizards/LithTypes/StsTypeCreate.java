
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.LithTypes;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;

public class StsTypeCreate extends StsWizardStep
{
    public StsLithTypesWizard wizard;
    public StsStatusPanel panel;
    private StsHeaderPanel header;
    private StsType type;

    public StsTypeCreate(StsWizard wizard)
    {
        super(wizard, new StsStatusPanel(), null, new StsHeaderPanel());

        this.wizard = (StsLithTypesWizard)wizard;
        panel = (StsStatusPanel) getContainer();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Type Definition");
        header.setSubtitle("Create Lithologic Type");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/LithLibrary.html");
    }

    public boolean start()
    {
        panel.setTitle("Creating Lithologic Type");
        run();
        return true;
    }

    public void run()
    {
        try
        {
            disableFinish();
            // Create New Type
            panel.setProgress(0.0f);
            panel.setText("Creating Lithologic Type");

//            type = new StsType(wizard.getTypeName(), wizard.getTypeColor(), wizard.getTypeTexture());

            panel.setText("Type Definition Complete");
            panel.setProgress(100.0f);

            enableFinish();
            disableCancel();
            actionManager.endCurrentAction();
            if(!success) return;
        }
        catch(Exception e)
        {
            success = false;
        }
    }

    public boolean end()
    {
        return true;
    }

}
