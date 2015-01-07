package com.Sts.Actions.Wizards.HandVelocity;

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
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
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

public class StsHandVelocitySelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew
{
    private StsFileTransferPanel selectionPanel;
    public HandVelocityFilenameFilter filenameFilter = null;

    String[] suffix = new String[] {"handvel", "HANDVEL", "dat"};  //allows for Focus VELDEF jobs SWC 9/8/09
    private DefaultListModel availableVolsListModel = new DefaultListModel();

    private JList availableObjectsList = new JList();

    private StsHandVelocityWizard wizard;
    private StsHandVelocitySelect wizardStep;

    private String currentDirectory = null;
    private StsAbstractFileSet availableFileSet;
    private StsAbstractFile[] availableFiles;

    StsModel model = null;
    private String outputDirectory = null;
    private StsFile currentSelectedFile = null;

    private StsButton viewButton = null;
    private StsComboBoxFieldBean vUnitsBean = null;
    private String velocityUnits = StsParameters.VEL_FT_PER_MSEC;
    private float scaleMultiplier = 1.0f;

    public StsHandVelocitySelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsHandVelocityWizard)wizard;
        this.wizardStep = (StsHandVelocitySelect)wizardStep;
        model = wizard.getModel();
        try
        {
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
			StsFilenameEndingFilter filenameFilter = new StsFilenameEndingFilter(suffix);
            selectionPanel = new StsFileTransferPanel( currentDirectory, filenameFilter, this, 350, 150);

            constructBeans();
            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public boolean hasDirectorySelection() { return true;  }
    public boolean hasReloadButton(){ return false;}
    public void    setReload(boolean reload) {}
    public boolean getReload() {return true;}
    public boolean hasArchiveItButton() { return false;}
    public void    setArchiveIt(boolean reload) {}
    public boolean getArchiveIt() { return true;}
    public boolean hasOverrideButton() { return false; }
    public void    setOverrideFilter(boolean override) {}
    public boolean getOverrideFilter() { return false; }

    public void fileSelected(StsAbstractFile selectedFile)
    {
        if(selectedFile == null)
            currentSelectedFile = null;
        if (selectedFile != currentSelectedFile)
            currentSelectedFile = (StsFile)selectedFile;
    }

    public void addFiles(StsAbstractFile[] files)
    {
        if( files.length <= 0) return;
        
        for (int n = 0; n < files.length; n++)
            wizard.addFile((StsFile)files[n]);
        fileSelected(files[0]);
    }
    
    public void availableFileSelected(StsAbstractFile selectedFile)
    {
    }
    
    public void addAllFiles(StsAbstractFile[] files)
    {
        for (int n = 0; n < files.length; n++)
            wizard.addFile((StsFile)files[n]);
        fileSelected(files[files.length-1]);
    }

    public void removeFiles(StsAbstractFile[] files)
    {
        for (int n = 0; n < files.length; n++)
            wizard.removeFile((StsFile)files[n]);
    }

    public void removeAllFiles()
    {
        wizard.removeFiles();
    }

    public void initialize()
    {
        filenameFilter = new HandVelocityFilenameFilter(suffix);
        setAvailableFiles();
        availableObjectsList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    }

    private void setAvailableFiles()
    {
        availableFileSet = StsFileSet.constructor(currentDirectory, filenameFilter);
        availableFiles = availableFileSet.getFiles();
        String[] fileEndings = getFilenameEndings( availableFiles);
        int cnt = 0;
        availableVolsListModel.removeAllElements();
        if(fileEndings != null)
        {
          for(int i = 0; i < fileEndings.length; i++)
          {
              if( model.getObjectWithName(StsSeismicVolume.class, fileEndings[i]) == null)
              {
                  availableVolsListModel.addElement(fileEndings[i]);
                  availableFiles[cnt] = availableFiles[i];
                  cnt++;
              }
          }
        }
    }

    private void constructBeans()
    {
        viewButton = new StsButton("View Selected File", "Review selected file contents", this, "showFile");
        vUnitsBean = new StsComboBoxFieldBean(this, "velocityUnits", "Velocity Units:", StsParameters.VEL_UNITS );
        vUnitsBean.setToolTipText("What units are the velocities in the selected file?");
	}

    public void showFile()
    {
        String line;
        StsAsciiFile selectedFile = new StsAsciiFile(currentSelectedFile);
        if(selectedFile == null)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "No file selected in list, must be in selected (left) list");
            return;
        }
        if(!selectedFile.openReadWithErrorMessage()) return;

        StsFileViewDialog dialog = new StsFileViewDialog(wizard.frame,"Hand Velocity File View", false);
        dialog.setVisible(true);
        dialog.setViewTitle("File - " + currentSelectedFile.getFilename());
        try {
            while (true)
            {
                line = selectedFile.readLine();
                if (line == null)
                    break;
                dialog.appendLine(line);
            }
        }
        catch (Exception ex)
        {
            StsMessageFiles.infoMessage("Unable to view selected file.");
        }
        finally
        {
            selectedFile.close();
        }
    }

    void jbInit() throws Exception
    {
        this.setLayout(new GridBagLayout());

        Color panelColor = Color.black;
        Font  panelFont = new java.awt.Font("Dialog", 1, 11);

        this.gbc.gridwidth = 2;
        addEndRow(selectionPanel);
        addEndRow(viewButton);

        this.gbc.gridwidth = 1;
        this.gbc.fill = this.gbc.HORIZONTAL;
        this.gbc.insets = new Insets(0,20,0,5);
        addEndRow(vUnitsBean);
    }

    public void setVelocityUnits(String unitsString)
    {
        velocityUnits = unitsString;
        StsProject project = wizard.getModel().getProject();
        scaleMultiplier = project.calculateVelScaleMultiplier(velocityUnits);
    }

    public String getVelocityUnits() { return velocityUnits; }
    public float getScaleMultiplier() { return scaleMultiplier; }
    public String getCurrentDirectory() { return currentDirectory; }

    public String[] getFilenameEndings(StsAbstractFile[] files)
    {
        int nFiles = files.length;
        String[] fileEndings = new String[nFiles];
        for(int n = 0; n < nFiles; n++)
            fileEndings[n] = filenameFilter.getFilenameStem(files[n].getFilename());
        return fileEndings;
    }

    final class HandVelocityFilenameFilter implements FilenameFilter
     {
         String[] filter;
         int length;

         public HandVelocityFilenameFilter(String[] filter)
         {
             this.filter = new String[filter.length];
             for (int i = 0; i < filter.length; i++)
                 this.filter[i] = filter[i].toLowerCase();
         }

         public boolean accept(File dir, String name)
         {
             for (int i = 0; i < filter.length; i++)
             {
                 if (name.toLowerCase().endsWith(filter[i]))
                     return true;
             }
             return false;
        }

        public String getFilenameStem(String filename)
        {
            int filterStart = 1;
            for (int i = 0; i < filter.length; i++)
            {
                if (filename.toLowerCase().indexOf(filter[i]) > 0)
                {
                    filterStart = filename.toLowerCase().indexOf(filter[i]);
                    break;
                }
            }
            return filename.substring(0, filterStart - 1);
		}
    }
}
