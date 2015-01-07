package com.Sts.Actions.Wizards.Pointset;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;

public class StsCreatePointSetWizard extends StsWizard
{
    public StsSetupPointSet setupVolume;
    public StsDefinePointSet definePointSet;

    private StsWizardStep[] mySteps =
    {
        definePointSet = new StsDefinePointSet(this),
        setupVolume = new StsSetupPointSet(this)
    };

    public StsCreatePointSetWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Define PointSet");
        dialog.getContentPane().setSize(500, 500);
        if(!super.start()) return false;
        return true;
    }

    public boolean end()
    {
        return super.end();
    }

    public void next()
    {
        gotoNextStep();
    }

    public void previous()
    {
        gotoPreviousStep();
    }
}
