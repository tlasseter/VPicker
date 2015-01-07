package com.Sts.Actions.Wizards.SeismicAttribute;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

import java.util.*;

public class StsSeismicAttributeWizard extends StsWizard implements Runnable
{
    private StsSelectVolume selectVolume;
    private StsDefineType defineType;
    private StsProcessCube processCube;

    private StsWizardStep[] mySteps =
    {
        defineType = new StsDefineType(this),
        selectVolume = new StsSelectVolume(this),
        processCube = new StsProcessCube(this)
    };

    public StsSeismicAttributeWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Seismic Attribute Calculation");
        dialog.getContentPane().setSize(500, 300);
        StsBoxSetSubVolumeClass boxSetClass = (StsBoxSetSubVolumeClass)model.getStsClass(StsBoxSetSubVolume.class);
        boxSetClass.setIsVisible(true);
        if(!super.start()) return false;
        disableFinish();
        return true;
    }

    public boolean end()
    {
        StsBoxSetSubVolumeClass boxSetClass = (StsBoxSetSubVolumeClass)model.getStsClass(StsBoxSetSubVolume.class);
        boxSetClass.setIsVisible(false);
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
//        if(currentStep == selectVolume) enableFinish();
        gotoNextStep();
    }

	public void cancel()
	{
        StsSeismicVolume volume = processCube.getVolume();
        if(volume != null) processCube.cancel();
		super.cancel();
	}


    public StsSeismicVolume getSeismicVolume()
    {
        return selectVolume.getVolume();
    }

    public ArrayList getAttributes() { return defineType.getAttributes(); }

	public boolean isDataFloat() { return selectVolume.isDataFloat(); }

    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsSeismicAttributeWizard msiVolumeWizard = new StsSeismicAttributeWizard(actionManager);
        msiVolumeWizard.start();
    }
}
