package com.Sts.Actions.Wizards.EclipseLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;

public class StsEclipseLoadWizard extends StsWizard
{
    public StsEclipseFileSelect fileSelect = null;
    public StsEclipseFileLoad fileLoad = null;

    private StsAbstractFile restartFile = null;
    private StsEclipseModel eclipseData;

    private StsWizardStep[] mySteps =
    {
        fileSelect = new StsEclipseFileSelect(this),
        fileLoad = new StsEclipseFileLoad(this)
    };

    public StsEclipseLoadWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Eclipse File Load");
        disableFinish();
        if(!super.start()) return false;
        return true;
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
        if(currentStep == fileSelect)
        {
            initializeDatasets();
            gotoLoad();
        }
        else
            gotoNextStep();
    }

    private void initializeDatasets()
    {
        restartFile = fileSelect.panel.getSelectedFile();
    }

    public void gotoLoad()
    {
        fileLoad.constructPanel();
        gotoStep(fileLoad);
    }

    public void initialize()
    {

    }

    public void addFile(StsFile file)
    {
        restartFile = file;
    }

    public void removeAllFiles() { restartFile = null; }

    public void removeFile(StsAbstractFile file)
    {
    	if(file == null) return;
        restartFile = null;
    }

    public StsAbstractFile getRestartFile() { return restartFile; }

    public StsEclipseModel getEclipseData()
    {
        return eclipseData;
    }

    public void setEclipseData(StsEclipseModel eclipseData)
    {
        this.eclipseData = eclipseData;
    }
}