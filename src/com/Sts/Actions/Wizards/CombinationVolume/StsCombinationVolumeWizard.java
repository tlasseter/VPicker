package com.Sts.Actions.Wizards.CombinationVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.SubVolume.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

import javax.swing.*;

public class StsCombinationVolumeWizard extends StsWizard
{
    final byte BOX_SET = 0;
    final byte SV_CUMMULATIVE = 1;
    final byte CUMMULATIVE = 2;
    final byte DIFFERENCE = 3;
    final byte MATH = 4;
    static final public String[] mvTypes = new String[] {"Box Set Summation", "Cumulative SubVolume", "Multi-PostStack3d Difference", "Individual Difference", "Multi-PostStack3d Math"};

    public StsDefineCombinationVolume defineVolume;
    public StsSetupCombinationVolume setupVolume;
    public StsSelectBoxSet selectBoxSet;
    public StsDefineBoxSetSubVolume defineBoxSet;
    public StsDefineCummulativeVolume cummVolume;
    public StsDefineCummulativeDifferenceVolume cummDiffVolume;
    public StsDefineDifferenceVolume diffVolume;
//    public StsDefinePointSet definePointSet;
    public StsDefineMathVolume mathVolume;

    private StsWizardStep[] mySteps =
    {
        defineVolume = new StsDefineCombinationVolume(this),
        selectBoxSet = new StsSelectBoxSet(this),
        defineBoxSet = new StsDefineBoxSetSubVolume(this),
        cummVolume = new StsDefineCummulativeVolume(this),
        cummDiffVolume = new StsDefineCummulativeDifferenceVolume(this),
        diffVolume = new StsDefineDifferenceVolume(this),
        mathVolume = new StsDefineMathVolume(this),
//        definePointSet = new StsDefinePointSet(this),
        setupVolume = new StsSetupCombinationVolume(this)
    };
    private StsSubVolume subVolume = null;

    public StsCombinationVolumeWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Define Combination PostStack3d");
        dialog.getContentPane().setSize(500, 500);
        if(!super.start()) return false;
        return true;
    }

    public boolean end()
    {
        return super.end();
    }

    public void previous()
    {
//        if(currentStep == definePointSet)
//            gotoStep(mathVolume);
        if(currentStep != defineVolume)
            gotoStep(defineVolume);
    }

    public void next()
    {
        if (currentStep == defineVolume)
        {
            switch(getType())
            {
                case BOX_SET:
                    gotoStep(selectBoxSet);
                    disableFinish();
                    break;
                case SV_CUMMULATIVE:
                    new StsMessage(model.win3d, StsMessage.WARNING, "Cumulative Sub-Volumes are not available yet.");
//                    gotoStep(cummVolume);
                    break;
                case CUMMULATIVE:
                    gotoStep(cummDiffVolume);
                    disableFinish();
                    break;
                case DIFFERENCE:
                    gotoStep(diffVolume);
                    disableFinish();
                    break;
                case MATH:
                    gotoStep(mathVolume);
                    disableFinish();
                    break;
            }
        }
        else if(currentStep == defineBoxSet)
        {
            gotoStep(setupVolume);
            enableFinish();
        }
        else if(currentStep == cummVolume)
        {
            gotoStep(setupVolume);
            enableFinish();
        }
        else if(currentStep == cummDiffVolume)
        {
            if(cummDiffVolume.getVolumes().length%2 != 0)
            {
                new StsMessage(this.frame, StsMessage.WARNING,"Must select an even number of volumes.");
                return;
            }
            StsSeismicVolume[] list = cummDiffVolume.getVolumes();
            for(int i=1; i<list.length; i++)
            {
                if(!list[0].sameAs(list[i]))
                {
                    new StsMessage(this.frame, StsMessage.WARNING,"PostStack3d " + (i+1) + " is not the same size as PostStack3d 1.");
                    return;
                }
            }
            gotoStep(setupVolume);
            enableFinish();
        }
        else if(currentStep == diffVolume)
        {

            if(!diffVolume.getVolumeOne().sameAs(diffVolume.getVolumeTwo()))
            {
                new StsMessage(this.frame, StsMessage.WARNING,"PostStack3d 1 is not the same size as PostStack3d 2.");
                return;
            }
            gotoStep(setupVolume);
            enableFinish();
        }
        else if(currentStep == mathVolume)
        {
            StsSeismicVolume[] list = mathVolume.getVolumes();
            for(int i=1; i<list.length; i++)
            {
                if(!list[0].sameAs(list[i]))
                {
                    new StsMessage(this.frame, StsMessage.WARNING,"PostStack3d " + (i+1) + " is not the same size as PostStack3d 1.");
                    return;
                }
            }
            gotoStep(setupVolume);
            enableFinish();
        }
        else
            gotoNextStep();
    }

    public boolean isInclusive() { return defineVolume.isInclusive(); }
    public byte getType() { return defineVolume.getType(); }
    public String getVolumeName() { return defineVolume.getVolumeName(); }
    public StsBoxSetSubVolume getBoxSet() { return (StsBoxSetSubVolume)selectBoxSet.getBoxSet(); }
    public void setSubVolume(StsSubVolume subVolume) { this.subVolume = subVolume; }
    public StsSubVolume getSubVolume() { return subVolume; }
    public JPanel getHistogramPanel() { return mathVolume.panel.getHistogramPanel(); }
}
