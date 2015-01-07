package com.Sts.Actions.Wizards.SimulationFile;

import com.Sts.Actions.Export.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

// import CFunction;
// import CPointer;

public class StsInitializeExportSimulationFile extends StsWizardStep implements Runnable
{
    StsStatusPanel status;
    StsHeaderPanel header;

    public int format = 0;
    public StsPropertyVolumeOld porosity = null;
    public StsPropertyVolumeOld permeability = null;
    public double kx, ky, kz;
    public String stemname;

    static public final int ECLIPSE = 1;
    static public final int RESCUE = 2;

    public StsInitializeExportSimulationFile(StsWizard wizard)
    {
    	super(wizard, new StsStatusPanel(), null, new StsHeaderPanel());
        status = (StsStatusPanel) getContainer();
        stemname = model.getName();
//        status.setLogFile(model.getLogFile());

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Simulation File Export");
        header.setSubtitle("Export File");
    }

	public boolean start()
    {
        status.setTitle("Creating the simulation file:");
		return true;
    }

    public void run()
    {
		success = createFile(status);
     //   wizard.enableFinish();
		if(success)
            status.setText("Simulation file created successfully.");
        else
        	status.setText("Unable to create simulation file.");

        status.sleep(1000);
//		wizard.finish();
    }


    public boolean end()
    {
		return success;
    }

    public void setFormat(int format) { this.format = format; }
    public void setStemname(String stemname) { this.stemname = stemname; }


	private boolean createFile(StsStatusPanel status)
    {
        this.status = status;
        boolean success = false;
        StsCursor cursor = new StsCursor(model.win3d, Cursor.WAIT_CURSOR);
        if(format == ECLIPSE) success = createEclipseFile(stemname, status);
        else if(format == RESCUE)
        {
            StsExportRescueFile file = new StsExportRescueFile(model, stemname, status);
            success = file.saveModel(status);
        }
        cursor.restoreCursor();
        return success;
    }


    private boolean createEclipseFile(String stemname, StsStatusPanel status)
    {
        try
        {
            status.setMaximum(100);
            float progress = 0;

            StsEclipseOutput eclipse = StsEclipseOutput.getInstance(model);
            eclipse.status = status;
            if(eclipse == null) return false;
            status.setText("Creating Eclipse file...");
            if(!eclipse.createEclipseFiles(stemname)) return false;
            return true;
        }
        catch(Exception e)
        {
            status.setProgress(0);
            StsException.outputException("StsInitializeExportSimulationFile.createFile() failed.", e, StsException.WARNING);
            return false;
        }
    }
}
