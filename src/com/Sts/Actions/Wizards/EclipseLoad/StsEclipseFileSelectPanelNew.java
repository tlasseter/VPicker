package com.Sts.Actions.Wizards.EclipseLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.DataTransfer.*;
import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsEclipseFileSelectPanelNew extends StsJPanel implements StsObjectTransferListener
{
	private StsEclipseLoadWizard wizard;
	private StsEclipseFileSelect wizardStep;
	private StsModel model;
	private String outputDirectory = null;
    private String currentDirectory = null;
    private StsAbstractFileSet fileset;
    private StsAbstractFileSetTransferPanel selectionPanel;
    private StsComboBoxFieldBean eclipseDataBean;

    static final String[] filterStrings = new String[]	{"unrst", "UNRST"};
    static final String jarFilename = "EclipseRestarts.jar";

    public StsEclipseFileSelectPanelNew(StsEclipseLoadWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = wizard;
		this.wizardStep = (StsEclipseFileSelect) wizardStep;
		this.model = wizard.getModel();
		try
		{
			StsModel model = wizard.getModel();
			if (model != null)
			{
				StsProject project = model.getProject();
				currentDirectory = project.getRootDirString();
				outputDirectory = project.getDataFullDirString();
			}
			else
			{
				currentDirectory = System.getProperty("user.dirNo"); // standalone testing
				outputDirectory = currentDirectory;
			}
            fileset = initializeFileSet();
			selectionPanel = new StsAbstractFileSetTransferPanel("Restart load", currentDirectory, this, 300, 100, fileset);
            eclipseDataBean = new StsComboBoxFieldBean(wizard, "eclipseData", "S2S-Eclipse Dataset", model.getObjectList(StsEclipseModel.class));
            constructPanel();

		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

    void constructPanel() throws Exception
	{
		this.gbc.fill = GridBagConstraints.BOTH;
		this.gbc.anchor = GridBagConstraints.NORTH;
		this.gbc.weighty = 1.0;
		add(selectionPanel);
        this.gbc.fill = GridBagConstraints.NONE;
        add(eclipseDataBean);
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
            StsFilenameEndingFilter filter = new StsFilenameEndingFilter(filterStrings);
            return fileset = StsFileSet.constructor(currentDirectory, filter);
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
        return true;
    }

    public StsAbstractFile getSelectedFile()
	{
        return (StsAbstractFile)(selectionPanel.getSelectedObject());
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
        return (StsAbstractFile[]) StsMath.trimArray(viewableFiles, nViewableFiles);
	}

    /** currently we can only view/edit seismic autopick file formats.  We need to be able at least view other
     *  ascii and binary files.
     * @param file being checked
     * @return true if viewable
     */
    private boolean isFileViewable(StsAbstractFile file)
    {
        return false;
    }
 
    public int getSelectedSurfaceIndex()
    {
        return selectionPanel.getSelectedIndices()[0];
    }

    public String getCurrentDirectory() { return currentDirectory; }
}
