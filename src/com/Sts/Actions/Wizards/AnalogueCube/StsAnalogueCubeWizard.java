package com.Sts.Actions.Wizards.AnalogueCube;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Reflect.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

public class StsAnalogueCubeWizard extends StsWizard implements Runnable
{
    private StsSelectVolume selectVolume;
    private StsDefineScope defineScope;
    private StsDefineType defineType;
    private StsProcessCube processCube;

    private StsAnalogueVolumeConstructor analogueVolumeConstructor = null;
    private StsSeismicVolume targetVolume = null;
	private String attributeVolumeName;
    private float correlThreshold = -1.0f;
    private byte optimizationMethod = SPIRAL_AVERAGE;
    private boolean hasExecuted = true;

    private Thread processThread = null;
    private boolean isDryRun = false;
	private StsProcessCubePanel processCubePanel = null;
    private StsProgressPanel progressPanel = null;

    private StsWizardStep[] mySteps =
    {
        defineType = new StsDefineType(this),
        selectVolume = new StsSelectVolume(this),
        defineScope = new StsDefineScope(this),
        processCube = new StsProcessCube(this)
    };

    public static final byte REAL = 0;
    public static final byte COMPLEX = 1;

    public static byte BYTE = 0;
    public static byte NATIVE = 1;

    public static final byte RUNNING_AVERAGE = 0;
    public static final byte SPIRAL_AVERAGE = 1;
    public static final byte SPIRAL_MAXIMUM = 2;

    public StsAnalogueCubeWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Aspect Energy Analogue Cube Analysis");
        dialog.getContentPane().setSize(500, 300);
		StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
		subVolumeClass.setIsVisible(false);
		subVolumeClass.setIsApplied(false);
        if(!super.start()) return false;
        return true;
    }

    public boolean end()
    {
		StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
		subVolumeClass.setIsVisible(false);
		subVolumeClass.setIsApplied(false);

        defineScope.clearBox();
        if(!hasExecuted) deleteLoadedAnalogueVolume();
       return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == selectVolume)
        {
            targetVolume = selectVolume.getVolume();
			String optionalName = selectVolume.getOptionalName();
			String stemname;
			if (optionalName.length() > 0)
				stemname = targetVolume.stemname + "." + optionalName + "." + StsAnalogueVolumeConstructor.ANALOGUE;
			else
				stemname = targetVolume.stemname + "." + StsAnalogueVolumeConstructor.ANALOGUE;

            StsSeismicVolume existingVolume = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, stemname);
            if (existingVolume != null) existingVolume.delete();
        }
        else if(currentStep == defineScope)
        {
            enableFinish();
        }
        gotoNextStep();
    }

/*
    public void dryRun()
    {
        Runnable dryRunRunnable = new Runnable()
        {
            public void run()
            {
                dryRun(processCube.panel);
            }
        };

        processThread = new Thread(dryRunRunnable);
        processThread.start();
    }
*/
    public void dryRun(StsProcessCubePanel panel, boolean isDryRun)
    {
        hasExecuted = false;
		this.processCubePanel = panel;
        progressPanel = processCubePanel.progressPanel;
		this.isDryRun = isDryRun;

        if(!isDryRun)
        {
            if(analogueVolumeConstructor != null)
            {
                endDataPlaneCalculation();
            }
            processCubePanel.buttonsReinitialize();
            return;
        }

		Runnable dryRunRunnable = new Runnable()
		{
			public void run()
			{
				runDryRun();
				processCubePanel.buttonsReinitialize();
			}
		};

		Thread dryRunThread = new Thread(dryRunRunnable);
		dryRunThread.start();
    }

	public void runDryRun()
	{
        initializeVolume(processCubePanel, false, true);
        if(analogueVolumeConstructor == null)
        {
            new StsMessage(null, StsMessage.WARNING, "Failed to construct/classInitialize analogue volume");
            return;
        }
        processCubePanel.displayHeader("Preview");
        analogueVolumeConstructor.initializeRun();

        startDataPlaneCalculation();
        byte[] byteArray = new byte[0];
        Class byteArrayClass = byteArray.getClass();
        StsMethod dryRunMethod = new StsMethod(analogueVolumeConstructor, "dryRun", new Class[] { byteArrayClass, Integer.class, Integer.class });
        StsSeismicVolume attributeVolume = getAttributeVolume();
        attributeVolume.setCalculateDataMethod(dryRunMethod);
		try
		{
			attributeVolume.initialize(model, "rw");
		}
		catch(Exception e)
		{
			StsMessage.printMessage("Failed to classInitialize analogueCubeVolume. Error: " + e.getMessage());
			attributeVolume.delete();
			return;
		}

        model.setCurrentObjectDisplayAndToolbar(attributeVolume);
//        model.setCurrentObject(attributeVolume);
//        model.win3d.getCursor3d().clearTextureDisplays();
        model.win3dDisplay();
    }

    public void fullRun(StsProcessCubePanel panel)
    {
		this.processCubePanel = panel;
        Runnable fullRunRunnable = new Runnable()
        {
            public void run()
            {
                fullRun();
            }
        };

        processThread = new Thread(fullRunRunnable);
        processThread.start();
    }

    public void fullRun()
    {
        progressPanel = processCubePanel.progressPanel;
        initializeVolume(processCubePanel, true, false);
		analogueVolumeConstructor.useByteCubes(false);
        if(analogueVolumeConstructor == null)
        {
            new StsMessage(null, StsMessage.WARNING, "Failed to construct/classInitialize analogue volume");
            return;
        }
        endDataPlaneCalculation();

        hasExecuted = analogueVolumeConstructor.fullRun();
        if(!hasExecuted)
            processCubePanel.buttonsReinitialize();
        else
		{
            processCubePanel.buttonsReinitialize();
			progressPanel.finished();
			progressPanel.appendLine("Initializing display of " + analogueVolumeConstructor.getVolume().getName());
		}
        StsSeismicVolume attributeVolume = getAttributeVolume();
		if(attributeVolume == null) return;
		if (!attributeVolume.initialize(model))
		{
			attributeVolume.delete();
			return;
		}
        model.setCurrentObject(attributeVolume);
        model.viewObjectChanged(this, attributeVolume);
        model.win3dDisplay();
        progressPanel.appendLine("Processing finished...");
    }


    public StsSubVolume getSubVolume()
    {
        return selectVolume.getSubVolume();
    }

    public void startDataPlaneCalculation()
    {
        StsSeismicVolume volume = getAttributeVolume();
        if (volume != null) volume.setCheckWriteCursorPlaneFlag(true);
    }

    public void endDataPlaneCalculation()
    {
        StsSeismicVolume volume = getAttributeVolume();
        if (volume != null) volume.setCheckWriteCursorPlaneFlag(false);
    }

    private void initializeVolume(StsProcessCubePanel panel, boolean createFiles, boolean useByteCubes)
    {
        StsSeismicVolume targetVolume = selectVolume.getVolume();
        StsSeismicVolume sourceVolume = defineScope.getSourceVolume();
        StsSubVolume targetSubVolume = selectVolume.getSubVolume();
		String optionalName = selectVolume.getOptionalName();
        boolean isDataFloat = selectVolume.isDataFloat();
        byte type = defineType.getType();
        byte spiralMethod = processCube.panel.getOptimizationMethod();

        StsBoxSubVolume sourceSubVolume = defineScope.getSourceBox();

//        deleteLoadedAnalogueVolume(targetVolume);
//        analogueVolumeConstructor = StsAnalogueVolumeConstructor.constructor(model, sourceVolume, targetVolume,
//            sourceSubVolume, targetSubVolume, resolution, type, spiralMethod, panel, createFiles);
        if(analogueVolumeConstructor == null)
            analogueVolumeConstructor = StsAnalogueVolumeConstructor.constructor(model, sourceVolume, targetVolume,
                sourceSubVolume, targetSubVolume, isDataFloat, type, spiralMethod, progressPanel, false, optionalName);
        if(createFiles && !analogueVolumeConstructor.checkCreateFiles())
		{
			analogueVolumeConstructor = null;
			return;
		}

		analogueVolumeConstructor.useByteCubes(useByteCubes);
		analogueVolumeConstructor.fillSourceCube();
        panel.displayHeader("Full Run");
        if (analogueVolumeConstructor == null)
        {
            new StsMessage(null, StsMessage.WARNING, "Failed to construct/classInitialize analogue volume");
            return;
        }
        getCorrelThreshold();
        analogueVolumeConstructor.initializeRun();
    }

    public float getCorrelThreshold()
    {
        correlThreshold = processCube.getCorrelationThreshold();
        if(analogueVolumeConstructor != null)
            analogueVolumeConstructor.setCorrelationThreshold(correlThreshold);
        return correlThreshold;
    }

    public void setCorrelThreshold(float threshold)
    {
        this.correlThreshold = threshold;
        if(analogueVolumeConstructor == null) return;

        analogueVolumeConstructor.setCorrelationThreshold(correlThreshold);
        StsSeismicVolume attributeVolume = analogueVolumeConstructor.getVolume();
        attributeVolume.deleteExistingFiles();
        attributeVolume.close();
        processCube.panel.buttonsReinitialize();
        processCube.panel.resetDryRunBtn();
    }

    public void deleteLoadedAnalogueVolume()
    {
		if(analogueVolumeConstructor == null) return;
        StsSeismicVolume attributeVolume = analogueVolumeConstructor.getVolume();
        if (attributeVolume != null)  attributeVolume.delete();
    }

    public StsAnalogueVolumeConstructor getAnalogueVolume() { return analogueVolumeConstructor; }

    public StsSeismicVolume getTargetVolume() { return targetVolume; }

    public StsSeismicVolume getAttributeVolume()
    {
        if(analogueVolumeConstructor != null) return analogueVolumeConstructor.getVolume();
        else return null;
    }

	public void killRun()
	{
		if(analogueVolumeConstructor != null) analogueVolumeConstructor.canceled = true;
	}
    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsAnalogueCubeWizard msiVolumeWizard = new StsAnalogueCubeWizard(actionManager);
        msiVolumeWizard.start();
    }
}
