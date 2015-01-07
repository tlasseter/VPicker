package com.Sts.Actions.Wizards.FracSetGolderLoad;

import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.IO.StsAbstractFile;
import com.Sts.IO.StsFile;
import com.Sts.IO.StsFilenameEndingFilter;
import com.Sts.Interfaces.StsFileTransferObjectFaceNew;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.StsProject;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.StsFileTransferPanel;

import java.awt.*;


/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsGolderFracSetSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew
{
	private StsFileTransferPanel selectionPanel;

	private StsGolderFracSetLoadWizard wizard;
	private StsGolderFracSetSelect wizardStep;
	private StsModel model;
	private String outputDirectory = null;
    private boolean archiveIt = true;


	public StsGolderFracSetSelectPanel(StsGolderFracSetLoadWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = wizard;
		this.wizardStep = (StsGolderFracSetSelect) wizardStep;
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
			String[] filterStrings = new String[]	{"csv", "txt", "Csv", "Txt", "CSV", "TXT", "FAB", "fab"};
			StsFilenameEndingFilter filter = new StsFilenameEndingFilter(filterStrings);

			selectionPanel = new StsFileTransferPanel(currentDirectory, filter, this, 300, 100, false);
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
		this.setLayout(new GridBagLayout());

		this.gbc.fill = GridBagConstraints.BOTH;
		this.gbc.anchor = GridBagConstraints.NORTH;
		this.gbc.weighty = 1.0;
		add(selectionPanel);
		this.gbc.fill = GridBagConstraints.HORIZONTAL;
		this.gbc.weighty = 0.0;
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
        StsFile[] files = wizard.getFracSetFiles();
        for(int i=0; i<files.length; i++)
            wizard.removeFile(files[i]);
    }

    public void fileSelected(StsAbstractFile selectedFile) { ; }
    public void removeFiles(StsAbstractFile[] files)
    {
        for(int i=0; i<files.length; i++)
           wizard.removeFile((StsFile)files[i]);
	}

	public boolean hasReloadButton()  	{ return false;	}
	public void setReload(boolean reload) {}
	public boolean getReload() 	{ return true;	}
    public boolean hasArchiveItButton() { return false; }
    public void setArchiveIt(boolean reload)  { archiveIt = reload; }
    public boolean getArchiveIt() { return archiveIt; }
    public boolean hasOverrideButton() { return false; }
    public void setOverrideFilter(boolean override) {}
    public boolean getOverrideFilter() { return false; }
    
}
