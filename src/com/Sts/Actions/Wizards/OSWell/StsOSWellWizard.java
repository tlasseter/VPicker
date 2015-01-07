//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.OSWell;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.OpenSpirit.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsOSWellWizard extends StsWizard
{
    int nWells;
    StsPoint[] topPoints = null;
    public StsOpenSpiritImport ospImport = null;
    public StsOSWellDatastore currentDatastore = null;

    static final int NUM_FILE_TYPES = 5;
    public static final byte HEADER = 0;
    public static final byte SURVEY = 1;
    public static final byte LOGS = 2;
    public static final byte TD = 3;
    public static final byte TOPS = 4;

    StsOSWellSelect selectWells = null;
    StsOSWellLoad loadWells = null;
    StsWizardStep[] mySteps;

    byte vUnits = StsParameters.DIST_NONE;
    byte hUnits = StsParameters.DIST_NONE;
    byte binaryHorzUnits = StsParameters.DIST_NONE;
    byte binaryVertUnits = StsParameters.DIST_NONE;

    static public final boolean debug = false;

    public StsOSWellWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 700);

        System.setSecurityManager(null);

        frame.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        ospImport = new StsOpenSpiritImport(this.model);
        if(!ospImport.initializeOpenSpirit())
        {
            frame.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            return;
        }
        // see if we have a persisted OpenSpirit well datastore and if not, make a new one.
        currentDatastore = (StsOSWellDatastore)model.getCreateStsClass(StsOSWellDatastore.class).getCurrentObject();
        if (currentDatastore == null)
        {
        	// create an empty datastore instance
        	currentDatastore = new StsOSWellDatastore(ospImport);
        }
        else
            currentDatastore.setOpenSpiritImport(ospImport);
        frame.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
        addSteps(true);
    }

    private void addSteps(boolean displayReloadAscii)
         {
             selectWells = new StsOSWellSelect(this);
             loadWells = new StsOSWellLoad(this);
             mySteps = new StsWizardStep[] {selectWells, loadWells};
             addSteps(mySteps);
         }

    public boolean checkPrerequisites()
              {
                  String ospHome = com.openspirit.OpenSpiritFactory.getOpenSpiritHome().toString();
                  if (ospHome == null)
                  {
                      reasonForFailure = "The OSP_HOME Environment Variable has not been set. Are you a licensed OpenSpirit user?\n\n";
                      return false;
                  }
                  return true;
              }

    public boolean start()
                   {

                       System.runFinalization();
                       System.gc();
                       dialog.setTitle("Load OpenSpirit Wells");
                       initialize();
                       this.disableFinish();
                       return super.start();
                   }

    public void initialize()
                        {
                            hUnits = model.getProject().getXyUnits();
                            vUnits = model.getProject().getDepthUnits();
                        }

    public StsOSWell[] getSelectedWells()
                             {
                                 return selectWells.getSelectedWells();
                             }

    public boolean end()
                                  {
                                      disConnectOpenSpirit();

                                      if (success)
                                      {
                                          model.setActionStatus(getClass().getName(), StsModel.STARTED);
                                      }
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

    public void finish()
                             {
                                 success = true;
                                 super.finish();
                             }

    public void cancel()
                        {
                            disConnectOpenSpirit();
                            if(currentDatastore != null)
                            {
                                currentDatastore.delete();
                                currentDatastore = null;
                                ospImport = null;
                            }
                            super.cancel();
                        }
    /**
     * Get the StsOpenSpiritImport instance.
                    * @return the StsOpenSpiritImport instance.
                    */
                   public StsOpenSpiritImport getOpenSpiritImport()
                   {
                       return ospImport;
                   }

    /**
     * Get the StsOSWellDatastore instance.
               * @return the StsOSWellDatastore instance.
               */
              public StsOSWellDatastore getWellDatastore()
              {
                  return currentDatastore;
              }

    /**
     * Disconnect from the OpenSpirit instance.
          */
         public void disConnectOpenSpirit()
         {
             ospImport.exitOpenSpirit();
         }

    static public void main(String[] args)
    {
        StsModel model = StsModel.constructor();
        StsActionManager actionManager = new StsActionManager(model);
        StsOSWellWizard wellWizard = new StsOSWellWizard(actionManager);
        wellWizard.start();
    }

}
