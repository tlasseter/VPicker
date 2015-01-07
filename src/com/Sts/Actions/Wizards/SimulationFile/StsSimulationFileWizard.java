package com.Sts.Actions.Wizards.SimulationFile;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Export.StsEclipseOutput;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;


public class StsSimulationFileWizard extends StsWizard
{
	private static final int SELECT_FILE_AND_FORMAT = 0;
	private static final int SELECT_PROPERTY_VOLUME = 1;
	private static final int EDIT_KXYZ = 2;
	private static final int EXPORT_FILE = 3;

	private StsSelectSimulationFile selectFile = new StsSelectSimulationFile(this);
	private StsSelectPropertyVolume selectProperty = new StsSelectPropertyVolume(this);
//	private StsEditKxyz editKxyz = new StsEditKxyz(this);
	private StsInitializeExportSimulationFile initializeFile = new StsInitializeExportSimulationFile(this);
    StsDefineMinCellSize cellSize = new StsDefineMinCellSize(this);
    StsOutputExportSimulationFile outputFile = new StsOutputExportSimulationFile(this);
    private StsWizardStep[] mySteps = { selectFile, selectProperty, initializeFile, cellSize, outputFile};

    String filename;
    int format;
    StsPropertyVolumeOld porosity, permeability;
    double kx, ky, kz;

	public StsSimulationFileWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
        disableFinish();
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Export Simulation File");
    	return super.start();
    }

    public boolean checkStartAction()
    {
        StsEclipseModel eclipseModel = (StsEclipseModel)model.getCurrentObject(StsEclipseModel.class);
        if(eclipseModel == null) return true;
        boolean yes = StsYesNoDialog.questionValue(model.win3d,"Eclipse Model already built.\n  Do you wish to delete it and rebuild?.");
        if(yes) eclipseModel.deleteModel(model);
        return yes;
    }


    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == selectFile)
        {
            filename = selectFile.getFilename();
            if(filename == null) return;
            initializeFile.setStemname(filename);
            format = selectFile.getFormat();
            initializeFile.setFormat(format);
            if(hasProperties()) gotoNextStep();
            else gotoStep(initializeFile);
        }
        else if(currentStep == selectProperty)
        {
            porosity = selectProperty.getSelectedPorosity();
            permeability = selectProperty.getSelectedPermeability();
            kx = selectProperty.getKx();
            ky = selectProperty.getKy();
            kz = selectProperty.getKz();
            initializeFile.porosity = porosity;
            initializeFile.permeability = permeability;
            initializeFile.kx = kx;
            initializeFile.ky = ky;
            initializeFile.kz = kz;
            gotoNextStep();
        }
        else
            gotoNextStep();
    }

    public void finish()
    {
    	next();
        super.finish();
    }


    public boolean end()
    {
        model.win3d.win3dDisplay();
        System.runFinalization();
        System.gc();
		return super.end();
    }

    private boolean hasProperties()
    {
        StsPropertyVolumeOld[] properties = (StsPropertyVolumeOld[])model.getCastObjectList(StsPropertyVolumeOld.class);
        if (properties == null || properties.length < 1 )
        {
            statusArea.setText ("No properties available: will build grid only.");
            return false;
        }
        else
            return true;
    }

    public void cancel()
    {
        StsEclipseOutput.clearTransients();
    }
}
