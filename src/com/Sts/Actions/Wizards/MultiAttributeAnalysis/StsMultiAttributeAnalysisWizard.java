package com.Sts.Actions.Wizards.MultiAttributeAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

public class StsMultiAttributeAnalysisWizard extends StsWizard
{
    private StsAttributeSelection attributeSelection;

    private StsWizardStep[] mySteps =
    {
        attributeSelection = new StsAttributeSelection(this)
    };

    public StsMultiAttributeAnalysisWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Define Multi-Attribute Vector");
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
        return "MAVName";
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

}
