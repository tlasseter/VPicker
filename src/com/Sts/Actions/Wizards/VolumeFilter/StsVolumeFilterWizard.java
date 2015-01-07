package com.Sts.Actions.Wizards.VolumeFilter;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

public class StsVolumeFilterWizard extends StsWizard implements Runnable
{
    private StsSelectVolume selectVolume;
    public StsDefineType defineType;
    private StsProcessCube processCube;
    
    public static final byte ANALYSIS_NONE = 0;
    public static final byte ANALYSIS_MEDIAN = 1;
    public static final byte ANALYSIS_MEAN = 2;
    

    private StsWizardStep[] mySteps =
    {
        selectVolume = new StsSelectVolume(this),
        defineType = new StsDefineType(this),       
        processCube = new StsProcessCube(this)
    };

    public StsVolumeFilterWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Seismic Volume Filter");
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

    //public ArrayList getAttributes() { return defineType.getAttributes(); }

	public boolean isDataFloat() { return selectVolume.isDataFloat(); }

    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsVolumeFilterWizard msiVolumeWizard = new StsVolumeFilterWizard(actionManager);
        msiVolumeWizard.start();
    }
}
