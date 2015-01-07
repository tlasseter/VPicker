package com.Sts.Actions.Wizards.FaultSticks;

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

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.WizardComponents.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.DataTransfer.*;
import com.Sts.DBTypes.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsFaultSticksSelectPanel extends StsJPanel implements StsObjectTransferListener
{
    private StsModel model = null;
    private String currentDirectory = null;
    private StsAbstractFileSet fileset;
    private StsAbstractFileSetTransferPanel selectionPanel;
    public StsUnitsGroupBox unitsPanel;
    public StsVerticalUnitsGroupBox verticalUnitsPanel;
    //private JPanel selectedObjectInfoPanel = new JPanel();

    private StsFaultSticksWizard wizard;
    private StsFaultSticksSelect wizardStep;

    static final String jarFilename = "faultSticks.jar";

    public StsFaultSticksSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.model = wizard.getModel();
        this.wizard = (StsFaultSticksWizard)wizard;
        this.wizardStep = (StsFaultSticksSelect)wizardStep;

        try
        {
            constructPanel();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    void constructPanel() throws Exception
    {
        unitsPanel = new StsUnitsGroupBox(model);
        verticalUnitsPanel = new StsVerticalUnitsGroupBox(model);
        model = wizard.getModel();
         if (model == null)
             currentDirectory = System.getProperty("user.dirNo"); // standalone testing
         else
             currentDirectory = model.getProject().getRootDirString();
         fileset = initializeFileSet();
         selectionPanel = new StsAbstractFileSetTransferPanel("Fault Stick Files Selector", currentDirectory, this, 600, 200, fileset);
         add(selectionPanel);
         gbc.weighty = 0.0;
         add(unitsPanel);
         add(verticalUnitsPanel);
    }

    public void initialize()
    {
        StsModel model = wizard.getModel();
        if(model != null)
            currentDirectory = wizard.getModel().getProject().getRootDirString();
        else
            currentDirectory = System.getProperty("user.dirNo"); // standalone testing

        selectionPanel.setCurrentDirectory(currentDirectory);
    }

    private StsAbstractFileSet initializeFileSet()
     {
         // Get files from a webstart jar
         if(Main.isWebStart && Main.isJarFile)
              return StsWebStartJar.constructor(jarFilename);
          // Load from jar files
          else if(Main.isJarFile)
              return StsJar.constructor(currentDirectory, jarFilename);
          // Load from ASCII/Binary files
          else
         {
             String binaryDirectory = model.getProject().getBinaryFullDirString();
             return StsAsciiBinaryFileSet.constructor(currentDirectory, binaryDirectory, StsFaultStickSet.fileGroups);
         }
     }

    public String getCurrentDirectory() { return currentDirectory; }

    public StsAbstractFile[] getSelectedFiles()
	{
        Object[] selectedObjects = selectionPanel.getSelectedObjects();
        int nFiles = selectedObjects.length;
        if(nFiles == 0) return new StsAbstractFile[0];
        StsAbstractFile[] selectedFiles = new StsAbstractFile[nFiles];
        for(int n = 0; n < nFiles; n++)
            selectedFiles[n] = (StsAbstractFile)selectedObjects[n];
        return selectedFiles;
	}

    public void addObjects(Object[] objects)
    {
        System.out.println("add Objects:");
        printObjects(objects);
        wizard.enableNext(selectionPanel.getNSelectedObjects() > 0);
    }

    private void printObjects(Object[] objects)
    {
        for (int n = 0; n < objects.length; n++)
            System.out.println("    " + objects[n].toString());
    }

    public void removeObjects(Object[] objects)
    {
        System.out.println("remove Objects:");
        printObjects(objects);
        wizard.enableNext(selectionPanel.getNSelectedObjects() > 0);
    }

    public void objectSelected(Object selectedObject)
    {
        System.out.println("selected Object:" + selectedObject.toString());
    }

    public boolean addAvailableObjectOk(Object object)
    {
        StsAbstractFile file = (StsAbstractFile)object;
        String name = file.name;
        return model.getObjectWithName(StsFaultStickSet.class, name) == null;
    }
}
