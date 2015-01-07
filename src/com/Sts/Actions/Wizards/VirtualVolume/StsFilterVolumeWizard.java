package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

public class StsFilterVolumeWizard extends StsWizard
{
    private StsFilterVVolume filterVirtualVolume;

    private StsWizardStep[] mySteps =
    {
        filterVirtualVolume = new StsFilterVVolume(this)
    };

    public StsFilterVolumeWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);

        filterVirtualVolume.panel.setSeismicVolumes(getSeismicVolumes());
        filterVirtualVolume.panel.setVirtualVolumes(getVirtualVolumes());
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Define Virtual PostStack3d");
        dialog.getContentPane().setSize(400, 200);
        if(!super.start()) return false;

        checkAddToolbar();

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
        gotoNextStep();
    }
    public String getVolumeName()
    {
        return filterVirtualVolume.getName();
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