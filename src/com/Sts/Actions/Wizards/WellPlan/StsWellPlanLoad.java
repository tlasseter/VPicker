package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

import java.awt.*;

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
public class StsWellPlanLoad extends StsWizardStep // implements Runnable  commented out because run method does GUI actions: nono on non-event thread. TJL 1/28/07
{
	public StsProgressPanel panel;
	private StsHeaderPanel header;
	StsLoadWellPlanWizard wizard;
	StsFile[] selectedFiles;

	public StsWellPlanLoad(StsWizard wizard)
	{
		super(wizard);
		this.wizard = (StsLoadWellPlanWizard)wizard;
		panel = StsProgressPanel.constructor(5, 50);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("Well Plan Selection");
		header.setSubtitle("Load Well Plan(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#WellPlan");
		header.setInfoText(wizardDialog, "(1) Once complete, press the Finish Button to dismiss the screen.");
	}
	public boolean start()
	 {
		 run();
		 return true;
	 }

	 public void run()
	 {
		 try
		 {
			 panel.appendLine("Loading selected files...");

             disablePrevious();

			 // turn off redisplay
			 model.disableDisplay();

			 // turn on the wait cursor
             StsGLPanel3d glPanel3d = model.win3d.getGlPanel3d();
             StsCursor cursor = new StsCursor(glPanel3d, Cursor.WAIT_CURSOR);
			 completeWellLoad();
			 cursor.restoreCursor();
			 model.win3d.cursor3d.initialize();
			 model.win3d.cursor3dPanel.setSliderValues();
			 model.getGlPanel3d().setDefaultView();
  //           glPanel.getActionManager().endCurrentAction();
			 model.enableDisplay();
			 model.win3dDisplay(); // display the wells
		 }
		 catch (Exception e)
		 {
			 success = false;
		 }
	 }

	 public void setSelectedFiles(StsFile[] selectedFiles)
	 {
		 this.selectedFiles = selectedFiles;
	 }

	 /*public void completeWellLoad()
	 {
		 for(int n = 0; n < selectedFiles.length; n++)
		 {
			 String filename = selectedFiles[n].createFilename();
			 StsKeywordIO.parseAsciiFilename(filename);
			 StsWellPlanSet wellPlanSet = new StsWellPlanSet(StsKeywordIO.name);
			 StsWellPlan wellPlan = new StsWellPlan(false);
			 if(!StsDBFile.readDatabaseObjectFile(selectedFiles[n].getPathname(), wellPlan, null))
			 {
				 new StsMessage(model.win3d,  StsMessage.WARNING, "Failed to load file " + filename);
				 return;
			 }
//            StsWellPlan wellPlan = (StsWellPlan)StsToolkit.deserializeObject(selectedFiles[n]);
			 wellPlan.setIndex(-1);
			 wellPlan.addToModel();
			 wellPlan.isConstructing = false;
			 wellPlanSet.addWellPlan(wellPlan);
			 StsPlatform platform = wellPlan.getPlatform();
			 if(platform != null)
			 {
				 StsPlatform existingPlatform = (StsPlatform)model.getObjectWithName(StsPlatform.class, platform.getName());
				 if(existingPlatform == null)
				 {
					 platform.setIndex(-1);
					 platform.addToModel();
				 }
			 }
		 }
		 disableCancel();
		 enableFinish();
		 success = true;
	 }*/

	 public void completeWellLoad()
	 {
         panel.initialize(selectedFiles.length);
         for(int n = 0; n < selectedFiles.length; n++)
		 {
			 String filename = selectedFiles[n].getFilename();
			 StsKeywordIO.parseAsciiFilename(filename);
             StsWellPlanSet wellPlanSet = new StsWellPlanSet(StsKeywordIO.name);
             StsWellPlan wellPlan = (StsWellPlan)StsDBFileObjectTrader.importStsObject(selectedFiles[n].getPathname(), null);
             if(wellPlan == null)
			 {
				 new StsMessage(model.win3d,  StsMessage.WARNING, "Failed to load file " + filename);
				 continue;
			 }
 			wellPlan.setIndex(-1);
			wellPlan.addToModel();
            wellPlan.isConstructing = false;
            wellPlanSet.addWellPlan(wellPlan);
			StsPlatform platform = wellPlan.getPlatform();
			if(platform != null)
			{
				StsPlatform existingPlatform = (StsPlatform)model.getObjectWithName(StsPlatform.class, platform.getName());
				if(existingPlatform == null)
				{
					platform.setIndex(-1);
					platform.addToModel();
				}
			}
            panel.setValue(n+1);
            panel.setDescription("Loaded well plan " + n + " of " + selectedFiles.length);
         }
         panel.appendLine("All well plans loaded successfully.\n");
         panel.setDescription("Loading Complete");


         disableCancel();
         enableFinish();
		 success = true;
	 }

     public boolean end()
     {
         return success;
     }
}
