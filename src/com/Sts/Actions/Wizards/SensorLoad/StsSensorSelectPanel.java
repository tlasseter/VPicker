package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
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

public class StsSensorSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew
{
	private StsFileTransferPanel selectionPanel;

	private StsSensorLoadWizard wizard;
	private StsSensorSelect wizardStep;
	private StsModel model;
	private String outputDirectory = null;
    private boolean archiveIt = true;
    
    private byte sensorType = StsSensorClass.OTHER;
    private boolean computeAttributes = true;
    private boolean autoConfigure = true;

	private StsSensor currentSensor = null;
	
	private StsJPanel typePanel = new StsGroupBox();
	private StsComboBoxFieldBean typeBean;
    private StsBooleanFieldBean computeAttsBean;
    private StsBooleanFieldBean autoConfigureBean;

	public StsSensorSelectPanel(StsSensorLoadWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = wizard;
		this.wizardStep = (StsSensorSelect) wizardStep;
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
			String[] filterStrings = new String[]	{"csv", "txt", "src"};
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
        typeBean = new StsComboBoxFieldBean(this, "typeString", "Sensor Type:", StsSensorClass.typeStrings);
        computeAttsBean = new StsBooleanFieldBean(this, "computeAttributes", "Compute Attributes?");
        computeAttsBean.setToolTipText("Will generate additional attributes on load, including cumulative magnitude and event count.");
        autoConfigureBean = new StsBooleanFieldBean(this, "autoConfigure", "Automatically Configure?");
        autoConfigureBean.setToolTipText("Will automatically determine number of stages and set colors when true.");
		typePanel.gbc.gridwidth = 2;
        typePanel.addEndRow(typeBean);
        typePanel.gbc.gridwidth = 1;
        typePanel.addToRow(computeAttsBean);
        typePanel.addEndRow(autoConfigureBean);
		add(typePanel);
	}

    public boolean hasDirectorySelection() { return true;  }

    public void addFiles(StsAbstractFile[] files)
    {
        for(int i=0; i<files.length; i++)
        {
           if(!wizard.addFile((StsFile)files[i]))
           {
               selectionPanel.removeFile(files[i]);
           }

        }
    }
    
    public void availableFileSelected(StsAbstractFile selectedFile)
    {
    }
    
    public void removeAllFiles()
    {
        StsSensorFile[] files = wizard.getSensorFiles();
        for(int i=0; i<files.length; i++)
            wizard.removeFile(files[i].file);
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
    
    public void setTypeString(String typeStg)
    {
    	sensorType = StsSensorClass.getTypeForName(typeStg);
    }
    public String getTypeString() { return StsSensorClass.typeStrings[sensorType]; }
    public byte getSensorType() { return sensorType; }
    public boolean getAutoConfigure() { return autoConfigure; }
    public void setAutoConfigure(boolean val)
    {
        autoConfigure = val;
    }
    public boolean getComputeAttributes() { return computeAttributes; }
    public void setComputeAttributes(boolean val)
    {
        computeAttributes = val;
    }
}
