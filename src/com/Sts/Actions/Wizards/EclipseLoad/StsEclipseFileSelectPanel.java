package com.Sts.Actions.Wizards.EclipseLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsEclipseFileSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew
{
	private StsFileTransferPanel selectionPanel;

	private StsEclipseLoadWizard wizard;
	private StsEclipseFileSelect wizardStep;
	private StsModel model;
	private String outputDirectory = null;
    private boolean archiveIt = true;

	public StsEclipseFileSelectPanel(StsEclipseLoadWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = wizard;
		this.wizardStep = (StsEclipseFileSelect) wizardStep;
		this.model = wizard.getModel();
		try
		{
			String currentDirectory;
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
			String[] filterStrings = new String[]	{"unrst", "UNRST"};
			StsFilenameEndingFilter filter = new StsFilenameEndingFilter(filterStrings);

			selectionPanel = new StsFileTransferPanel(currentDirectory, filter, this, 300, 100, true);
            jbInit();

		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

    public void initialize()
    {
        wizard.rebuild();
    }

    void jbInit() throws Exception
	{
		this.gbc.fill = GridBagConstraints.BOTH;
		this.gbc.anchor = GridBagConstraints.NORTH;
		this.gbc.weighty = 1.0;
		add(selectionPanel);
	}

    public boolean hasDirectorySelection() { return true;  }

    public void addFiles(StsAbstractFile[] files)
    {
        for(int i=0; i<files.length; i++)
           wizard.addFile((StsFile)files[i]);
    }

    public void availableFileSelected(StsAbstractFile selectedFile)
    {
    }

    public void removeAllFiles()
    {
        wizard.removeAllFiles();
    }

    public void fileSelected(StsAbstractFile selectedFile) { ; }
    public void removeFiles(StsAbstractFile[] files)
    {
        for(int i=0; i<files.length; i++)
           wizard.removeFile(files[i]);
	}

	public boolean hasReloadButton()  	{ return false;	}
	public void setReload(boolean reload) {}
	public boolean getReload() 	{ return true;	}
    public boolean hasArchiveItButton() { return true; }
    public void setArchiveIt(boolean reload)  { archiveIt = reload; }
    public boolean getArchiveIt() { return archiveIt; }
    public boolean hasOverrideButton() { return false; }
    public void setOverrideFilter(boolean override) {}
    public boolean getOverrideFilter() { return false; }
}