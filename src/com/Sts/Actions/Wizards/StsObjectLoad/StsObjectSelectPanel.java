package com.Sts.Actions.Wizards.StsObjectLoad;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import java.awt.*;
import java.io.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsObjectSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew
{
	
	private StsFileTransferPanel selectionPanel;

	private StsObjectLoadWizard wizard;
	private StsObjectSelect wizardStep;
	private StsModel model;
	private String outputDirectory = null;

	public StsObjectSelectPanel(StsObjectLoadWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = wizard;
		this.wizardStep = (StsObjectSelect) wizardStep;
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

			StsFilenameFilterObj filter = new StsFilenameFilterObj("obj");
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
        StsFile[] files = wizard.getObjectFiles();
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
    public void setArchiveIt(boolean reload)  { }
    public boolean getArchiveIt() { return false; }
    public boolean hasOverrideButton() { return false; }
    public void setOverrideFilter(boolean override) {}
    public boolean getOverrideFilter() { return false; }
    
    class StsFilenameFilterObj extends StsFilenameFilter
    {
    	public StsFilenameFilterObj(String fmt)
        {
            this.format = fmt;
        }

    	public boolean accept(File dir, String filename)
    	{
    		StsWellKeywordIO.parseBinaryFilename(filename);
    		return StsWellKeywordIO.format.equalsIgnoreCase(format);
    	}
    	
    	public String getType(String filename)
    	{
    		StsWellKeywordIO.parseBinaryFilename(filename);
    		return StsWellKeywordIO.group;
    	}
    }
}
