package com.Sts.Actions.Wizards.MultiAttributeAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

public class StsMultiAttrAnalysisWizard extends StsWizard
{
    private StsAttributeSelection attributeSelection;

    private StsWizardStep[] mySteps =
    {
        attributeSelection = new StsAttributeSelection(this)
    };

    public StsMultiAttrAnalysisWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }
    
    public boolean checkPrerequisites()
    {
    	int nSeismicVolumes = ((StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class)).length;
        StsVirtualVolumeClass vvClass = (StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class);
    	int nVirtualVolumes = vvClass.getVirtualVolumes().length;
    	if((nSeismicVolumes + nVirtualVolumes) < 2) 
    	{
    		reasonForFailure = "Must have at least two volumes to run multi-attribute analysis. Volumes can be virtual or real. Please load more volumes or run the Virtual Volume Step from the Seismic Attributes Workflow first.";
    		return false;
    	}
    	else
    		return true;
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
