package com.Sts.Actions.Wizards.SubVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

public class StsSubVolumeWizard extends StsWizard
{
    private StsDefineSubVolume defineSubVolume;
    private StsSetupSubVolume setupSubVolume;
    private StsDefineSurfaceSubVolume surfaceSubVolume;
    private StsDefineDualSurfaceSubVolume dualSurfaceSubVolume;
    private StsSelectBoxSet selectBoxSetSubVolume;
    private StsDefineBoxSetSubVolume defineBoxSetSubVolume;
    private byte offsetDomain;

    boolean isCanceled = false;

    private StsWizardStep[] mySteps =
    {
        defineSubVolume = new StsDefineSubVolume(this),
        surfaceSubVolume = new StsDefineSurfaceSubVolume(this),
        dualSurfaceSubVolume = new StsDefineDualSurfaceSubVolume(this),
        selectBoxSetSubVolume = new StsSelectBoxSet(this),
        defineBoxSetSubVolume = new StsDefineBoxSetSubVolume(this),
        setupSubVolume = new StsSetupSubVolume(this)
    };

    private StsSubVolume subVolume = null;

    public StsSubVolumeWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Define Sub-Volume");
        dialog.getContentPane().setSize(500, 500);
        if(!super.start()) return false;

        return true;
    }

    public boolean end()
    {
        return super.end();
    }

    public void cancel()
    {
        if(subVolume != null && !subVolume.delete()) return;
        isCanceled = true;
        super.cancel();
    }

    public boolean isCanceled()
    {
        return isCanceled;
    }
    public void previous()
    {
        if(currentStep != defineSubVolume)
        {
            gotoStep(defineSubVolume);
        }
    }

    public void next()
    {
        if (currentStep == defineSubVolume)
        {
            switch(getSubVolumeType())
            {
                case StsSubVolume.BOX_SET:
                    gotoStep(selectBoxSetSubVolume);
                    break;
                case StsSubVolume.SINGLE_SURFACE:
                    if(model.getNObjects(StsModelSurface.class) < 1)
                        new StsMessage(frame, StsMessage.WARNING, "Need at least one horizon to define this sub-volume type.\n" +
                        		"If surfaces have already been loaded, please run the Define->Horizons Wizard in the Velocity Modeling Workflow.");
                    else
                    {
                        gotoStep(surfaceSubVolume);
                        enableFinish();
                        disableNext();
                    }
                    break;
                case StsSubVolume.DUAL_SURFACE:
                    if(model.getNObjects(StsModelSurface.class) < 2)
                        new StsMessage(frame, StsMessage.WARNING, "Need at least two horizons to define this sub-volume type.\n" +
                		"If surfaces have already been loaded, please run the Define->Horizons Wizard in the Velocity Modeling Workflow.");
                    else
                    {
                        gotoStep(dualSurfaceSubVolume);
                        enableFinish();
                        disableNext();
                    }
                    break;
                case StsSubVolume.WELL_SET:
//                    if(model.getNObjects(StsWell.class) < 1)
//                        new StsMessage(frame, StsMessage.WARNING, "Need at least one well to define this sub-volume type.");
//                    else
                        new StsMessage(frame, StsMessage.WARNING, "Well constrained sub-volumes are not currently supported.");
                    break;
                case StsSubVolume.RESERVOIR_UNIT:
                        new StsMessage(frame, StsMessage.WARNING, "Reservoir unit constrained sub-volumes are not currently supported.");
                    break;

            }

        }
        else if (currentStep == selectBoxSetSubVolume)
        {
            StsBoxSetSubVolume boxSetSubVolume = (StsBoxSetSubVolume)subVolume;
            defineBoxSetSubVolume.setBoxSet(boxSetSubVolume);
            gotoStep(defineBoxSetSubVolume);
            enableFinish();
            disableNext();
        }
    }

    public void finish()
    {
        model.subVolumeChanged();
        super.finish();
        return;
    }

    public void setSubVolume(StsSubVolume subVolume) { this.subVolume = subVolume; }
    public StsSubVolume getSubVolume() { return subVolume; }

    public byte getSubVolumeType()
    {
        return defineSubVolume.getSubVolumeType();
    }
    
    public byte getOffsetDomain()
    {
        return offsetDomain;
    }

    public void setOffsetDomain(byte offsetDomain)
    {
        this.offsetDomain = offsetDomain;
    }

    public String getSubVolumeName() { return defineSubVolume.getSubVolumeName(); }
    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsSubVolumeWizard subVolumeWizard = new StsSubVolumeWizard(actionManager);
        subVolumeWizard.start();
    }
}
