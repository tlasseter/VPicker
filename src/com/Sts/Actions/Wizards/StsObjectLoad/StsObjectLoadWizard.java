
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.StsObjectLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

public class StsObjectLoadWizard extends StsWizard // implements Runnable
{
    public StsObjectSelect objectSelect = null;
    public StsObjectLoad objectLoad = null;
    public StsFile[] files = null;
    
    private StsWizardStep[] mySteps =
    {
        objectSelect = new StsObjectSelect(this),
        objectLoad = new StsObjectLoad(this)
        
    };

    public StsObjectLoadWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 600); 
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("StsObject Load Wizard");
        initialize();
        disableFinish();
        if(!super.start()) return false;
        return true;
    }

    public void initialize()
    {
    	;
    }

    public boolean end()
    {
        model.setActionStatus(getClass().getName(), StsModel.STARTED);    	
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

    public void addFile(StsFile objFile)
    {
        files = (StsFile[]) StsMath.arrayAddElement(files, objFile);
    }

    public void removeFile(StsFile objFile)
    {
    	if(objFile != null)
    		files = (StsFile[])StsMath.arrayDeleteElement(files, objFile);
    }

    public StsFile[] getObjectFiles() { return files; }
}
