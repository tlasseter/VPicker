package com.Sts.Actions.Wizards.VolumeInterpolation;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

public class StsVolumeInterpolationWizard extends StsWizard
{
    private StsDefinition definitionStep;

    private StsWizardStep[] mySteps =
    {
    		definitionStep = new StsDefinition(this)
    };

    public StsVolumeInterpolationWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Interpolate Volume");
        dialog.getContentPane().setSize(400, 200);
        if(!super.start()) return false;

        checkAddToolbar();

        return true;
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        gotoNextStep();
    }

    public String getVolumeName()
    {
    	String name = "InterpolatedVolume";
    	if(getVolume() != null)
    		name = getVolume().getName() + "_Interpolated";
        return name;
    }

    public StsSeismicVolume[] getSeismicVolumes()
    {
        return (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
    }

    public StsVirtualVolume[] getVirtualVolumes()
    {
        StsVirtualVolumeClass vvClass = (StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class);
        return vvClass.getVirtualVolumes();
    }

    public StsSeismicVolume getVolume()
    {
    	return (StsSeismicVolume)definitionStep.getSelectedVolume();
    }
}
