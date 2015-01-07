package com.Sts.Actions.Wizards.PreStackExport;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;

public class StsPreStackExportWizard extends StsWizard
{
    private StsExportType exportType;
//    private StsSelectData selectData;
//    private StsExportData exportData;
    
    public static byte TWO_D = 0;
    public static byte THREE_D = 1;
    
    public static byte VELOCITY = 0;
    public static byte STACK = 1;
    public static byte SEMBLANCE = 2;
    public static byte GATHER_DATA = 3;
    
    private StsWizardStep[] mySteps =
    {
        exportType = new StsExportType(this),
//        selectData = new StsSelectData(this),
//        exportData = new StsExportData(this)
    };

    public StsPreStackExportWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Export Pre-Stack Data");
        dialog.getContentPane().setSize(400, 300);
        if(!super.start()) return false;
        
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
    
    public int getExportType()
    {
        return exportType.getExportType();
    }
    public int getExportDimension()
    {
        return exportType.getExportDimension();
    }
}
