package com.Sts.Actions.Wizards.FracSetMVLoad;

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

public class StsMVFracSetSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 2705137639441940182L;

	private StsFileTransferPanel selectionPanel;

	private StsMVFracSetLoadWizard wizard;
	private StsMVFracSetSelect wizardStep;
	private StsModel model;
	private String outputDirectory = null;
    private boolean archiveIt = true;
    
	
	public StsMVFracSetSelectPanel(StsMVFracSetLoadWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = wizard;
		this.wizardStep = (StsMVFracSetSelect) wizardStep;
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
			String[] filterStrings = new String[]	{"csv", "txt", "Csv", "Txt", "CSV", "TXT"};
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
        StsMVFracSetFile[] files = wizard.getFracSetFiles();
        for(int i=0; i<files.length; i++)
            wizard.removeFile(files[i].file);
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
