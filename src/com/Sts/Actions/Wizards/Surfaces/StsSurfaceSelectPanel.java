package com.Sts.Actions.Wizards.Surfaces;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.WizardComponents.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.DataTransfer.*;
import com.Sts.Utilities.*;
import com.Sts.DBTypes.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSurfaceSelectPanel extends StsJPanel implements StsObjectTransferListener
{
    private StsSurfaceWizard wizard;
    private StsSurfaceSelect wizardStep;

    private StsModel model = null;
    private String currentDirectory = null;
    private StsAbstractFileSet fileset;
    private StsAbstractFileSetTransferPanel selectionPanel;
    boolean reloadAscii = false;
    // Use these when we have something to show in the surface info panel
    /*
    private StsGroupBox selectedSurfaceInfoPanel = new StsGroupBox();
    private JTextPane surfaceInfoLbl = new JTextPane();
    */
    private StsGroupBox parametersPanel = new StsGroupBox("Selection Parameters");
    private StsBooleanFieldBean reloadAsciiBean = new StsBooleanFieldBean(this, "reloadAscii", "Reload files from ascii");

    public StsUnitsGroupBox unitsPanel;
    public StsVerticalUnitsGroupBox verticalUnitsPanel;

    static final String jarFilename = "grids.jar";

    public StsSurfaceSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        model = wizard.getModel();
        this.wizard = (StsSurfaceWizard)wizard;
        this.wizardStep = (StsSurfaceSelect)wizardStep;

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
        parametersPanel.gbc.fill = gbc.HORIZONTAL;
        parametersPanel.add(reloadAsciiBean);
        gbc.fill = gbc.BOTH;
        gbc.weighty = 0.0;
        add(parametersPanel);
		if (model == null)
			currentDirectory = System.getProperty("user.dirNo"); // standalone testing
		else
			currentDirectory = model.getProject().getRootDirString();

        fileset = initializeFileSet();
		selectionPanel = new StsAbstractFileSetTransferPanel("Surface Selector", currentDirectory, this, 600, 200, fileset);

        gbc.weighty = 1.0;
        add(selectionPanel);

        gbc.weighty = 0.0;
        add(unitsPanel);
        add(verticalUnitsPanel);
        /*
        // Leave info panel off for now until we decide what info to put in it
        selectedSurfaceInfoPanel.gbc.fill = gbc.BOTH;
        surfaceInfoLbl.setText("Selected Surface Information");
		surfaceInfoLbl.setEditable(false);
		surfaceInfoLbl.setBackground(Color.lightGray);
		surfaceInfoLbl.setFont(new Font("Dialog", 0, 10));

        selectedSurfaceInfoPanel.setMinimumSize(new Dimension(600, 50));
        selectedSurfaceInfoPanel.addEndRow(surfaceInfoLbl);

        gbc.weighty = 1.0;
        add(selectedSurfaceInfoPanel);
        */
	}

    public void initialize()
    {
        StsModel model = wizard.getModel();
        if(model != null)
            currentDirectory = wizard.getModel().getProject().getRootDirString();
        else
            currentDirectory = System.getProperty("user.dirNo"); // standalone testing

        selectionPanel.setCurrentDirectory(currentDirectory);
        wizard.enableNext(false);
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
            return StsAsciiBinaryFileSet.constructor(currentDirectory, binaryDirectory, StsSurface.fileGroups);
        }
    }

    StsAbstractFile getCurrentFile()
    {
        StsAbstractFile currentFile = (StsAbstractFile)selectionPanel.getSelectedObject();
        if(currentFile != null) return currentFile;
        return (StsAbstractFile)selectionPanel.getSelectedObjects()[0];
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
        return model.getObjectWithName(StsSurface.class, name) == null;
    }

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

    public StsAbstractFile[] getViewableFiles()
	{
        Object[] selectedObjects = selectionPanel.getSelectedObjects();
        int nFiles = selectedObjects.length;
        if(nFiles == 0) return new StsAbstractFile[0];
        StsAbstractFile[] viewableFiles = new StsAbstractFile[nFiles];
        int nViewableFiles = 0;
        for(int n = 0; n < nFiles; n++)
        {
            StsAbstractFile selectedFile = (StsAbstractFile)selectedObjects[n];
            if(isFileViewable(selectedFile))
                viewableFiles[nViewableFiles++] = (StsAbstractFile)selectedObjects[n];
        }
        return (StsAbstractFile[])StsMath.trimArray(viewableFiles, nViewableFiles);
	}

    /** currently we can only view/edit seismic autopick file formats.  We need to be able at least view other
     *  ascii and binary files.
     * @param file being checked
     * @return true if viewable
     */
    private boolean isFileViewable(StsAbstractFile file)
    {
        String format = file.format;
        if(format == null) return false;
        if(!format.equals(StsAbstractFile.asciiFormat)) return false;
        return file.group.equals(StsSurface.seismicGrp);
    }

    public String[] getSelectedFilenames()
    {
        StsAbstractFile[] files = getSelectedFiles();
        int nFiles = files.length;
        if(nFiles == 0) return new String[] { "none" };
        String[] filenames = new String[nFiles];
        for(int n = 0; n < nFiles; n++)
            filenames[n] = files[n].filename;
        return filenames;
    }

    public String[] getSelectedSurfaceNames()
    {
        StsAbstractFile[] files = getSelectedFiles();
        int nFiles = files.length;
        if(nFiles == 0) return new String[] { "none" };
        String[] surfaceNames = new String[nFiles];
        for(int n = 0; n < nFiles; n++)
            surfaceNames[n] = files[n].name;
        return surfaceNames;
    }

    public Object[] getSelectedSurfaceObjects()
    {
        return getSelectedFiles();
    }

    public int getSelectedSurfaceIndex()
    {
        return selectionPanel.getSelectedIndices()[0];
    }

    public boolean getReloadAscii()
    {
        return reloadAscii;
    }
    public void setReloadAscii(boolean value)
    {
        if(reloadAscii == value) return;
        reloadAscii = value;
        ((StsAsciiBinaryFileSet)fileset).reloadAscii = value;
        selectionPanel.initializeSetAvailableObjects();
    }

    public String getCurrentDirectory() { return currentDirectory; }
}
