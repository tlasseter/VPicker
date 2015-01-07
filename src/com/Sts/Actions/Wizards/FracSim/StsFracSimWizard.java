package com.Sts.Actions.Wizards.FracSim;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

public class StsFracSimWizard extends StsWizard
{
    private StsDefineFracSimVolume defineFracSimVolume;
    private StsSetupFracSimVolume setupFracSimVolume;
    private StsDefineSingleSurfaceFracSimVolume singleSurfaceFracSimVolume;
    private StsDefineDualSurfaceFracSimVolume dualSurfaceFracSimVolume;

    private String fracSimType = SINGLE_SURFACE;
    private String name = "fracSim";

    boolean isCanceled = false;

    private StsWizardStep[] mySteps =
    {
        defineFracSimVolume = new StsDefineFracSimVolume(this),
        singleSurfaceFracSimVolume = new StsDefineSingleSurfaceFracSimVolume(this),
        dualSurfaceFracSimVolume = new StsDefineDualSurfaceFracSimVolume(this),
        setupFracSimVolume = new StsSetupFracSimVolume(this)
    };

    private StsFracSimVolume fracSimVolume = null;

	static final public String SINGLE_SURFACE = "Single Surface";
    static final public String DUAL_SURFACE = "Dual Surface";
    static final public String[] fracSimTypes = new String[] {SINGLE_SURFACE, DUAL_SURFACE };

    public StsFracSimWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Define Fracture Simulation Model");
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
        isCanceled = true;
        super.cancel();
    }

    public boolean isCanceled()
    {
        return isCanceled;
    }
    public void previous()
    {
        if(currentStep != defineFracSimVolume)
        {
            gotoStep(defineFracSimVolume);
        }
    }

    public void next()
    {
        if (currentStep == defineFracSimVolume)
        {
            if(fracSimType ==  SINGLE_SURFACE)
            {
                if(model.getNObjects(StsModelSurface.class) < 1)
                    new StsMessage(frame, StsMessage.WARNING, "Need at least one surface to define this fracture simulation model type.");
                else
                {
                    gotoStep(singleSurfaceFracSimVolume);
                    enableFinish();
                    disableNext();
                }
            }
            else if(fracSimType == DUAL_SURFACE)
            {
                if(model.getNObjects(StsModelSurface.class) < 2)
                    new StsMessage(frame, StsMessage.WARNING, "Need at least two surface to define this fracture simulation model.");
                else
                {
                    gotoStep(dualSurfaceFracSimVolume);
                    enableFinish();
                    disableNext();
                }
            }
        }
    }

    public void finish()
    {
        model.subVolumeChanged();
        super.finish();
        return;
    }

    public String getType() { return fracSimType; }
    public void setType(String type) { fracSimType = type; }

    public String getName() { return name; }
    public void setName(String name) { this.name = name; }

    public void setFracSimVolume(StsFracSimVolume fracSimVolume) { this.fracSimVolume = fracSimVolume; }
    public StsFracSimVolume getFracSimVolume() { return fracSimVolume; }

    static public void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsFracSimWizard fracSimWizard = new StsFracSimWizard(actionManager);
        fracSimWizard.start();
    }
}