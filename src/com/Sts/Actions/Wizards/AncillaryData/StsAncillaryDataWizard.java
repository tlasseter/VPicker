package com.Sts.Actions.Wizards.AncillaryData;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class StsAncillaryDataWizard extends StsWizard
{
    private StsFile[] ancillaryFiles = new StsFile[0];
    private StsAncillaryData[] ancillaryData = new StsAncillaryData[0];

    public StsAncillaryDataSelect ancillarySelect = new StsAncillaryDataSelect(this);
    public StsAncillaryDataDefine ancillaryDefine = new StsAncillaryDataDefine(this);
    public StsAncillaryAssign ancillaryAssign = new StsAncillaryAssign(this);
    public StsAncillaryDataLoad ancillaryLoad = new StsAncillaryDataLoad(this);

    private StsWizardStep[] mySteps =
        {
        ancillarySelect, ancillaryDefine, ancillaryAssign, ancillaryLoad
    };

    public StsAncillaryDataWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Load Ancillary Data (pdf, cgm, doc, txt,...)");
        dialog.getContentPane().setSize(500, 600);
        return super.start();
    }

    public boolean end()
    {
        if (success) model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == ancillarySelect)
        {
            createAncillaryData();
            gotoNextStep();
        }
        else if(currentStep == ancillaryDefine)
        {
            if(!ancillaryDefine.assignToWell())
            {
                ancillaryLoad.constructPanel();
                gotoStep(ancillaryLoad);
            }
            else
                gotoNextStep();
        }
        else if(currentStep == ancillaryAssign)
        {
            ancillaryLoad.constructPanel();
            gotoNextStep();
        }
    }

    public void finish()
    {
        super.finish();
    }

    public boolean createAncillaryData()
    {
       if(model == null)
           return false;

       if(ancillaryData.length != 0)
           ancillaryData = new StsAncillaryData[0];

       try
       {
           disablePrevious();

           int nSelected = ancillaryFiles.length;
           int nLoaded = 0;

           for(int n = 0; n < nSelected; n++)
           {
               if (ancillaryFiles[n] == null) continue;
               // See if it is already in the database.
               if(model.getObjectWithName(StsAncillaryData.class, ancillaryFiles[n].getFilename()) != null)
               {
                   boolean ans = StsYesNoDialog.questionValue(model.win3d,"File " + ancillaryFiles[n].getFilename() +
                                            " already in database.\n\nWould you like to replace it?");
                   if(ans)
                       model.getObjectWithName(StsAncillaryData.class, ancillaryFiles[n].getFilename()).delete();
                   else
                       continue;
               }
               makeAncillaryData(ancillaryFiles[nLoaded]);
               nLoaded++;
           }

           return true;
       }
       catch (Exception e)
       {
           StsException.outputException("StsAncillaryDataWizard.createAncillaryData() failed.", e, StsException.WARNING);
           return false;
       }
    }

    private boolean makeAncillaryData(StsAbstractFile file)
    {
        StsAncillaryData object = null;
        String filename = file.getPathname();
        byte type = StsAncillaryData.getTypeFromFilename(filename);
        object = new StsAncillaryData(filename, StsColor.GREEN, false, type);

        if(object != null)
            ancillaryData = (StsAncillaryData[]) StsMath.arrayAddElement(ancillaryData, object);
        return true;
    }

    public boolean loadAncillaryData(StsProgressPanel panel)
    {
        double progressValue = 0.0f;
        try
       {
           int nSelected = ancillaryData.length;
           int nLoaded = 0;
           panel.initialize(nSelected);
           for(int n = 0; n < nSelected; n++)
           {
               if(ancillarySelect.panel.getArchiveIt())
               {
                   if(!StsFile.copy(ancillaryData[n].getOriginalFilename(), ancillaryData[n].getFilename()))
                   {
                       panel.appendLine("\nSource file (" + ancillaryData[n].getOriginalFilename() + ") does not exist or cannot be copied to database.");
                       new StsMessage(model.win3d, StsMessage.ERROR, "Source file (" + ancillaryData[n].getOriginalFilename()
                                      + ") does not exist or cannot be copied to database.");
                       continue;
                   }
               }
               if (ancillaryData[n] == null)
                   continue;
               ancillaryData[nLoaded].addToModel();
               panel.appendLine("Successfully added " + ancillaryData[n].getName() + " to project.\n");
               nLoaded++;
               panel.setValue(n+1);
               panel.setDescription("Loaded " + nLoaded + " of " + nSelected);
           }
           panel.appendLine("Archival of ancillary data complete.");
           panel.setDescription("Load Complete");
           panel.finished();
           enableFinish();

           return true;
       }
       catch (Exception e)
       {
           panel.appendLine("StsAncillaryDataWizard.loadAncillaryData() failed.");
           panel.setDescriptionAndLevel("Exception thrown", StsProgressBar.ERROR);
           panel.finished();
           StsException.outputException("StsAncillaryDataWizard.loadAncillaryData() failed.", e, StsException.WARNING);
           return false;
       }
    }

    public void addFile(StsFile file, String outputDirectory)
    {
        ancillaryFiles = (StsFile[]) StsMath.arrayAddElement(ancillaryFiles, file);
    }

    public void removeFile(StsFile file)
    {
        ancillaryFiles = (StsFile[]) StsMath.arrayDeleteElement(ancillaryFiles, file);
    }

    public void removeFiles()
    {
        if(ancillaryFiles == null) return;
        for(int i = 0; i<ancillaryFiles.length; i++)
            ancillaryFiles = (StsFile[]) StsMath.arrayDeleteElement(ancillaryFiles, ancillaryFiles[i]);

        ancillaryFiles = null;
    }

    public StsFile[] getAncillaryDataFiles() { return ancillaryFiles; }
    public StsAncillaryData[] getAncillaryDataObjects() { return ancillaryData; }
}
