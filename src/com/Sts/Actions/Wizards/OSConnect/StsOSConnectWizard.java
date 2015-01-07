
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.OSConnect;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;

public class StsOSConnectWizard extends StsWizard
{
    public StsSelectServer selectServer = null;
    public StsServerConnect serverConnect = null;

    private StsWizardStep[] mySteps =
    {
        selectServer = new StsSelectServer(this),
        serverConnect = new StsServerConnect(this)
    };

    public StsOSConnectWizard(StsActionManager actionManager)
    {
        super(actionManager,500,500);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("OpenSpirit Connect");
        disableFinish();
        return super.start();
    }

    public boolean end()
    {
        return super.end();
    }

    public void previous()
    {
       gotoPreviousStep();
    }

    public void next()
    {
        gotoNextStep();
    }

    public void finish()
    {
        super.finish();
    }
    
    public boolean establishConnection(StsProgressTextPanel panel)
    {
    	return true;
    }
}
